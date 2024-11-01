use anyhow::{anyhow, Result};

use rust_htslib::bam;

use bio::io::fastq;

use bwa::BwaAligner;

use log::{debug, info, warn};

use std::collections::BTreeMap;

mod cli;

mod mapping_status;
use mapping_status::{get_mapping_status, MappingStatus};

mod read_pair_io;
use read_pair_io::{ReadPair, ReadPairIterator};

mod utils;
use utils::fastq_to_unmapped_fragments;
use utils::silence_stderr;

mod output;
use output::OutputWriter;

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

fn do_work(opts: &cli::CliOptions) -> Result<()> {
    info!(target: "IO", "Accessing Fastq read pairs from disk");
    let mut read_pair_iter =
        ReadPairIterator::new(opts.fastq_first.clone(), opts.fastq_second.clone())?;

    info!(target: "IO", "Creating output folders in {}", opts.output_folder.display());
    let mut output_writer = OutputWriter::new(opts.output_folder.clone())?;

    info!(target: "IO", "Loading BWA indices and creating aligners");
    let aligners = load_aligners(&opts.index)?;

    let mut total_alignments: usize = 0;
    loop {
        let mut read_pairs = read_pair_iter.batch(opts.batch_size);
        let n_read_pairs = read_pairs.len();

        if n_read_pairs == 0 {
            break;
        } else {
            info!(target: "BWA", "Aligning {} read pairs", read_pairs.len());
            for aligner in aligners.iter() {
                align_batch(aligner, &mut read_pairs, &opts)?;
            }
            for read_pair in read_pairs {
                use MappingStatus::*;
                match (read_pair.read1.status, read_pair.read2.status) {
                    (Mapped, Mapped) => continue,

                    (Suspicious, _) | (_, Suspicious) => {
                        warn!(target: "IO", "E7:Read pair {} has a suspicious alignment", read_pair.id);
                        continue;
                    }

                    (Unknown, _) | (_, Unknown) => {
                        warn!(target: "IO", "E8:Read pair {} has an unknown alignment status", read_pair.id);
                        continue;
                    }

                    (Unmapped, Unmapped) => {
                        output_writer.write("uuu1", &read_pair.read1.fastq)?;
                        output_writer.write("uuu2", &read_pair.read2.fastq)?;
                    }

                    (Mapped, Unmapped) => {
                        output_writer.write("mum1", &read_pair.read1.fastq)?;
                        output_writer.write("muu2", &read_pair.read2.fastq)?;
                    }

                    (Unmapped, Mapped) => {
                        output_writer.write("muu1", &read_pair.read1.fastq)?;
                        output_writer.write("mum2", &read_pair.read2.fastq)?;
                    }

                    (Mapped, PartiallyUnmapped(positions)) => {
                        output_writer.write("mpm1", &read_pair.read1.fastq)?;
                        output_writer.write("mpp2", &read_pair.read2.fastq)?;

                        let fragments =
                            fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions)?;
                        output_writer.write_batch("frag", &fragments)?;
                    }

                    (PartiallyUnmapped(positions), Mapped) => {
                        output_writer.write("mpp1", &read_pair.read1.fastq)?;
                        output_writer.write("mpm2", &read_pair.read2.fastq)?;
                        let fragments =
                            fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions)?;
                        output_writer.write_batch("frag", &fragments)?;
                    }

                    (Unmapped, PartiallyUnmapped(positions)) => {
                        output_writer.write("upu1", &read_pair.read1.fastq)?;
                        output_writer.write("upp2", &read_pair.read2.fastq)?;
                        let fragments =
                            fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions)?;
                        output_writer.write_batch("frag", &fragments)?;
                    }

                    (PartiallyUnmapped(positions), Unmapped) => {
                        output_writer.write("upp1", &read_pair.read1.fastq)?;
                        output_writer.write("upu2", &read_pair.read2.fastq)?;
                        let fragments =
                            fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions)?;
                        output_writer.write_batch("frag", &fragments)?;
                    }

                    (PartiallyUnmapped(positions1), PartiallyUnmapped(positions2)) => {
                        output_writer.write("ppp1", &read_pair.read1.fastq)?;
                        output_writer.write("ppp2", &read_pair.read2.fastq)?;
                        let fragments1 =
                            fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions1)?;
                        output_writer.write_batch("frag", &fragments1)?;
                        let fragments2 =
                            fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions2)?;
                        output_writer.write_batch("frag", &fragments2)?;
                    }
                }
            }
        }
        total_alignments += n_read_pairs;
        output_writer.flush()?;
        info!(target: "PROGRESS", "{} read pairs processed, {} reads written", total_alignments, output_writer.written());
    }
    info!(target: "SUMMARY", "{} read pairs processed, {} reads written, {} fragments written", total_alignments, output_writer.written(), output_writer.fragments());
    Ok(())
}

#[derive(Debug)]
struct MappedReadPair {
    #[allow(dead_code)]
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

fn align_batch(
    aligner: &BwaAligner,
    read_pairs: &mut [ReadPair],
    opts: &cli::CliOptions,
) -> Result<()> {
    let mut reads = Vec::new();
    for read_pair in read_pairs.iter() {
        match (&read_pair.read1.status, &read_pair.read2.status) {
            (MappingStatus::Mapped, MappingStatus::Mapped) => {
                continue;
            }
            (MappingStatus::Suspicious, _) | (_, MappingStatus::Suspicious) => {
                warn!("Read pair {} has a suspicious alignment", read_pair.id);
                continue;
            }
            _ => {
                reads.push(read_pair.read1.fastq.clone());
                reads.push(read_pair.read2.fastq.clone());
            }
        }
    }
    let alignments = align_reads(aligner, &reads, opts.threads.clone())?;

    for input_pair in read_pairs {
        let id = &input_pair.id;
        if alignments.contains_key(id) {
            let mapped_pair = alignments.get(id).unwrap();
            for record in mapped_pair.read1.iter() {
                input_pair.read1.alignments.push(record.clone());
            }
            for record in mapped_pair.read2.iter() {
                input_pair.read2.alignments.push(record.clone());
            }
            input_pair.read1.status = get_mapping_status(&input_pair.read1.alignments, &opts)?;
            input_pair.read2.status = get_mapping_status(&input_pair.read2.alignments, &opts)?;
        }
    }

    Ok(())
}

fn align_reads(
    aligner: &BwaAligner,
    reads: &[fastq::Record],
    threads: usize,
) -> Result<BTreeMap<String, MappedReadPair>> {
    debug!(target: "BWA", "Aligning {} reads", reads.len());
    let mut alignments = BTreeMap::new();
    let bwa_result = silence_stderr(|| aligner.align_fastq_records(&reads, true, true, threads))??;

    {
        // debugging
        use std::time::{SystemTime, UNIX_EPOCH};

        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards")
            .as_millis();

        let filename = format!("debug_{}.bam", timestamp);
        debug!("Writing BAM file: {}", filename);
        let header = aligner.create_bam_header();
        let mut writer = bam::Writer::from_path(&filename, &header, bam::Format::Bam)?;
        for aln in bwa_result.iter() {
            writer.write(aln)?;
        }
    }
    bwa_result.into_iter().for_each(|r| {
        let name = std::str::from_utf8(r.qname())
            .expect("Invalid UTF-8 in read name")
            .to_string();
        alignments
            .entry(name.clone())
            .or_insert_with(|| MappedReadPair::new(&name))
            .insert(r);
    });

    if alignments.len() == 0 {
        return Err(anyhow!("No alignments were generated"));
    }

    Ok(alignments)
}
