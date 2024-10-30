#![allow(unused_variables)]
#![allow(dead_code)]
use anyhow::{anyhow, Result};

use rust_htslib::bam;
use rust_htslib::bam::{record::Cigar, Format};

use bio::io::fastq;

use bwa::BwaAligner;
#[allow(unused_imports)]
use log::{debug, info, warn};

use std::collections::{HashMap, HashSet};

mod cli;

mod read_pair_io;
use read_pair_io::ReadPairIterator;

mod utils;
use utils::{block_filter, silence_stderr};

fn main() -> Result<()> {
    env_logger::init();

    let opts = cli::parse_cli()?;
    do_work(&opts)?;
    Ok(())
}
    
fn load_aligners(paths: &[std::path::PathBuf]) -> Result<Vec<BwaAligner>> {
    let mut aligners = Vec::<BwaAligner>::new();
    for pathbuf in paths.iter() {
        let aligner = silence_stderr(|| BwaAligner::from_path(pathbuf))??;
        aligners.push(aligner);
    }
    Ok(aligners)
}

#[derive(Debug)]
struct MappedReadPair {
    id: String,
    read1: Vec<bam::Record>,
    read2: Vec<bam::Record>,
}

impl MappedReadPair {
    fn new(id: &str) -> Self {
        MappedReadPair {
            id: id.to_string(),
            read1: Vec::new(),
            read2: Vec::new(),
        }
    }

    fn insert(&mut self, read: bam::Record) {
        if read.is_first_in_template() {
            self.read1.push(read);
        } else {
            self.read2.push(read);
        }
    }
}

fn align_reads(
    aligner: &BwaAligner,
    reads: &[fastq::Record],
    threads: usize,
) -> Result<HashMap<String, MappedReadPair>> {
    let mut alignments = HashMap::new();
    silence_stderr(|| aligner.align_fastq_records(&reads, true, true, threads))??
        .into_iter()
        .for_each(|r| {
            let name = std::str::from_utf8(r.qname())
                .expect("Invalid UTF-8 in read name")
                .to_string();
            alignments
                .entry(name.clone())
                .or_insert_with(|| MappedReadPair::new(&name))
                .insert(r);
        });
    Ok(alignments)
}

fn do_work(opts: &cli::CliOptions) -> Result<()> {
    info!(target: "IO", "Loading BWA indices and creating aligners");
    let aligners = load_aligners(&opts.index)?;
    let aligner = &aligners[0];
    let header = aligner.create_bam_header();

    info!(target: "IO", "Accessing Fastq read pairs from disk");
    let mut read_pair_iter =
        ReadPairIterator::new(opts.fastq_first.clone(), opts.fastq_second.clone())?;

    let mut writer = bam::Writer::from_stdout(&header, Format::Sam)?;
    let reads = read_pair_iter.batch(opts.batch_size);
    let mut count = 0;

    info!(target: "BWA", "Aligning {} reads", reads.len());
    if reads.len() > 0 {
        let alignments = align_reads(aligner, reads.as_slice(), opts.threads)?;
        for (id, mapped_pair) in alignments {
            for record in mapped_pair.read1.iter() {
                writer.write(record)?;
            }
            for record in mapped_pair.read2.iter() {
                writer.write(record)?;
            }
            {
                // debug only
                let id = &mapped_pair.id;
                let mapping_status1 = get_mapping_status(&mapped_pair.read1);
                let mapping_status2 = get_mapping_status(&mapped_pair.read2);
                if mapping_status1.is_err() || mapping_status2.is_err() {
                    return Err(anyhow!("Error processing unmapped positions"));
                }
                let mapping_status1 = mapping_status1.unwrap();
                let mapping_status2 = mapping_status2.unwrap();

                match (mapping_status1, mapping_status2) {
                    (MappingStatus::Mapped, MappingStatus::Mapped) => {
                        // Reject
                        count += 1;
                    }
                    (MappingStatus::Suspicious | MappingStatus::Unknown, _) | (_, MappingStatus::Suspicious | MappingStatus::Unknown) => {
                        // Error
                        warn!("Read pair {} has a suspicious alignment", id);
                    }
                    (MappingStatus::PartiallyUnmapped(_), _) | (_, MappingStatus::PartiallyUnmapped(_)) => {
                        // Allow
                        debug!("Read pair {} is partially unmapped", id);
                    }
                    _ => {
                        // Allow
                        debug!("Read pair {} is allowed", id);
                    }
                }
            }
        }
    } else {
        warn!(target: "IO", "No read pair found");
    }
    debug!("{} read pairs mapped in proper pairs", count);

    Ok(())
}

fn get_clipped_positions(record: &bam::Record) -> Vec<usize> {
    if record.is_unmapped() {
        return (0..record.seq_len()).collect();
    }

    let mut positions = Vec::new();
    let mut pos = 0; // keeps track of position in the read

    let cigar = record.cigar();
    let cigar_ops: Box<dyn Iterator<Item = &Cigar>> = if record.is_reverse() {
        Box::new(cigar.iter().rev())
    } else {
        Box::new(cigar.iter())
    };

    for cigar_op in cigar_ops {
        match cigar_op {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                positions.extend(pos..pos + *len as usize);
                pos += *len as usize;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                pos += *len as usize;
            }
            _ => {}
        }
    }
    positions
}

#[derive(Debug)]
enum MappingStatus {
    Mapped, // All positions map to the reference (indels + mismatches permitted, also split/supplementary reads are merged)
    Unmapped, // Read is designated unmapped by the mapper
    PartiallyUnmapped(HashSet<usize>), // Some positions map, while others are clipped
    Suspicious, // All positions are unmapped, but the mapper didn't designate it as unmapped
    Unknown, // Reads have not been aligned yet, or other ambiguous status
}

fn get_mapping_status(records: &[bam::Record]) -> Result<MappingStatus> {
    if records.is_empty() {
        return Err(anyhow!("No records to process"));
    }

    let read_name = records[0].qname();
    let first_or_second = records[0].is_first_in_template();
    let len = records[0].seq_len();

    if (records.len() == 1) && records[0].is_unmapped() {
        return Ok(MappingStatus::Unmapped);
    }

    if (records.len() > 1) && records[0].is_unmapped() {
        return Ok(MappingStatus::Suspicious);
    }

    let mut positions: HashSet<usize> =
        get_clipped_positions(&records[0]).iter().cloned().collect();

    for record in &records[1..] {
        if record.qname() != read_name {
            return Err(anyhow!("Read names do not match"));
        }

        if record.is_first_in_template() != first_or_second {
            return Err(anyhow!("Reads are not all the same position in the pair"));
        }

        if record.seq_len() != len {
            return Err(anyhow!("Reads are not all the same length"));
        }

        let p: HashSet<usize> = get_clipped_positions(record).iter().cloned().collect();

        if p.is_empty() {
            return Ok(MappingStatus::Mapped);
        }

        positions = positions.intersection(&p).cloned().collect();

        if positions.is_empty() {
            return Ok(MappingStatus::Mapped);
        }
    }
    if positions.is_empty() {
        return Ok(MappingStatus::Mapped);
    }

    if positions.len() == len {
        return Ok(MappingStatus::Suspicious);
    }

    // Mark read as partially unmapped if it passes block size and average block quality checks
    let positions = block_filter(positions, &records[0].qual(), 10, 20.0); // arbitrary values for min block size and quality

    if positions.is_empty() {
        return Ok(MappingStatus::Mapped);
    }

    Ok(MappingStatus::PartiallyUnmapped(positions))
}
