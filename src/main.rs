#![allow(unused_variables)]
#![allow(dead_code)]
use anyhow::Result;

use rust_htslib::bam;

use bio::io::fastq;

use bwa::BwaAligner;
#[allow(unused_imports)]
use log::{debug, info, warn};

use std::collections::BTreeMap;

mod cli;

mod mapping_status;
use mapping_status::{get_mapping_status, MappingStatus};

mod read_pair_io;
use read_pair_io::{ReadPair, ReadPairIterator};

mod utils;
use utils::silence_stderr;
use utils::fastq_to_unmapped_fragments;

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
    info!(target: "IO", "Loading BWA indices and creating aligners");
    let aligners = load_aligners(&opts.index)?;
    let header = aligners[0].create_bam_header();

    info!(target: "IO", "Accessing Fastq read pairs from disk");
    let mut read_pair_iter =
        ReadPairIterator::new(opts.fastq_first.clone(), opts.fastq_second.clone())?;

    info!(target: "IO", "Creating output folders in {}", opts.output_folder.display());
    std::fs::create_dir_all(&opts.output_folder)?;
    std::fs::create_dir_all(&opts.output_folder.join("unmapped_pairs"))?;
    std::fs::create_dir_all(&opts.output_folder.join("mapped_unmapped_pairs"))?;
    std::fs::create_dir_all(&opts.output_folder.join("mapped_partial_pairs"))?;
    std::fs::create_dir_all(&opts.output_folder.join("unmapped_partial_pairs"))?;
    std::fs::create_dir_all(&opts.output_folder.join("partial_pairs"))?;
    std::fs::create_dir_all(&opts.output_folder.join("unmapped_fragments"))?;

    // Open all output files
    let mut uuu1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_pairs/u1.fastq"))?;
    let mut uuu2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_pairs/u2.fastq"))?;
    let mut mum1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_unmapped_pairs/m1.fastq"))?;
    let mut muu2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_unmapped_pairs/u2.fastq"))?;
    let mut muu1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_unmapped_pairs/u1.fastq"))?;
    let mut mum2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_unmapped_pairs/m2.fastq"))?;
    let mut mpm1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_partial_pairs/m1.fastq"))?;
    let mut mpp2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_partial_pairs/p2.fastq"))?;
    let mut mpp1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_partial_pairs/p1.fastq"))?;
    let mut mpm2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("mapped_partial_pairs/m2.fastq"))?;
    let mut upu1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_partial_pairs/u1.fastq"))?;
    let mut upp2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_partial_pairs/p2.fastq"))?;
    let mut upp1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_partial_pairs/p1.fastq"))?;
    let mut upu2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_partial_pairs/u2.fastq"))?;
    let mut ppp1 = bio::io::fastq::Writer::to_file(opts.output_folder.join("partial_pairs/p1.fastq"))?;
    let mut ppp2 = bio::io::fastq::Writer::to_file(opts.output_folder.join("partial_pairs/p2.fastq"))?;
    let mut frag = bio::io::fastq::Writer::to_file(opts.output_folder.join("unmapped_fragments/fragments.fastq"))?;
    
    loop {
        let mut read_pairs = read_pair_iter.batch(opts.batch_size);

        if read_pairs.len() == 0 {
            break;
        } else {
            info!(target: "BWA", "Aligning {} read pairs", read_pairs.len());
            for aligner in aligners.iter() {
                align_batch(aligner, &mut read_pairs, opts.threads)?;
            }
            for read_pair in read_pairs {
                use MappingStatus::*;
                match (read_pair.read1.status, read_pair.read2.status) {

                    (Mapped, Mapped) => continue,

                    (Suspicious, _) | (_, Suspicious) => {
                        warn!(target: "IO", "E7:Read pair {} has a suspicious alignment", read_pair.id);
                        continue;
                    },

                    (Unknown, _) | (_, Unknown) => {
                        warn!(target: "IO", "E8:Read pair {} has an unknown alignment status", read_pair.id);
                        continue;
                    },

                    (Unmapped, Unmapped) => {
                        uuu1.write_record(&read_pair.read1.fastq)?;
                        uuu2.write_record(&read_pair.read2.fastq)?;
                    },

                    (Mapped, Unmapped) => {
                        mum1.write_record(&read_pair.read1.fastq)?;
                        muu2.write_record(&read_pair.read2.fastq)?;
                    },

                    (Unmapped, Mapped) => {
                        muu1.write_record(&read_pair.read1.fastq)?;
                        mum2.write_record(&read_pair.read2.fastq)?;
                    },

                    (Mapped, PartiallyUnmapped(positions)) => {
                        mpm1.write_record(&read_pair.read1.fastq)?;
                        mpp2.write_record(&read_pair.read2.fastq)?;
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions)? {
                            frag.write_record(&fragment)?;
                        }
                    },

                    (PartiallyUnmapped(positions), Mapped) => {
                        mpp1.write_record(&read_pair.read1.fastq)?;
                        mpm2.write_record(&read_pair.read2.fastq)?;
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions)? {
                            frag.write_record(&fragment)?;
                        }
                    },

                    (Unmapped, PartiallyUnmapped(positions)) => {
                        upu1.write_record(&read_pair.read1.fastq)?;
                        upp2.write_record(&read_pair.read2.fastq)?;
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions)? {
                            frag.write_record(&fragment)?;
                        }
                    },

                    (PartiallyUnmapped(positions), Unmapped) => {
                        upp1.write_record(&read_pair.read1.fastq)?;
                        upu2.write_record(&read_pair.read2.fastq)?;
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions)? {
                            frag.write_record(&fragment)?;
                        }
                    },

                    (PartiallyUnmapped(positions1), PartiallyUnmapped(positions2)) => {
                        ppp1.write_record(&read_pair.read1.fastq)?;
                        ppp2.write_record(&read_pair.read2.fastq)?;
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read1.fastq, &positions1)? {
                            frag.write_record(&fragment)?;
                        }
                        for fragment in fastq_to_unmapped_fragments(&read_pair.read2.fastq, &positions2)? {
                            frag.write_record(&fragment)?;
                        }
                    },
                }
            }
        }
        uuu1.flush()?;
        uuu2.flush()?;
        mum1.flush()?;
        muu2.flush()?;
        muu1.flush()?;
        mum2.flush()?;
        mpm1.flush()?;
        mpp2.flush()?;
        mpp1.flush()?;
        mpm2.flush()?;
        upu1.flush()?;
        upp2.flush()?;
        upp1.flush()?;
        upu2.flush()?;
        ppp1.flush()?;
        ppp2.flush()?;
        frag.flush()?;
    }
    uuu1.flush()?;
    uuu2.flush()?;
    mum1.flush()?;
    muu2.flush()?;
    muu1.flush()?;
    mum2.flush()?;
    mpm1.flush()?;
    mpp2.flush()?;
    mpp1.flush()?;
    mpm2.flush()?;
    upu1.flush()?;
    upp2.flush()?;
    upp1.flush()?;
    upu2.flush()?;
    ppp1.flush()?;
    ppp2.flush()?;
    frag.flush()?;

    Ok(())
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

fn align_batch(
    aligner: &BwaAligner,
    read_pairs: &mut [ReadPair],
    threads: usize) -> Result<()> {
    let mut reads = Vec::new();
    for read_pair in read_pairs.iter() {
        match (&read_pair.read1.status, &read_pair.read2.status) {
            (MappingStatus::Mapped, MappingStatus::Mapped) => {
                continue;
            },
            (MappingStatus::Suspicious, _) | (_, MappingStatus::Suspicious) => {
                warn!("Read pair {} has a suspicious alignment", read_pair.id);
                continue;
            },
            _ => {
                reads.push(read_pair.read1.fastq.clone());
                reads.push(read_pair.read2.fastq.clone());
            }
        }
    }
    let alignments = align_reads(aligner, &reads, threads)?;

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
            input_pair.read1.status = get_mapping_status(&input_pair.read1.alignments)?;
            input_pair.read2.status = get_mapping_status(&input_pair.read2.alignments)?;
        }
    }
        
    Ok(())
}

fn align_reads(
    aligner: &BwaAligner,
    reads: &[fastq::Record],
    threads: usize,
) -> Result<BTreeMap<String, MappedReadPair>> {
    debug!("Aligning {} reads", reads.len());
    let mut alignments = BTreeMap::new();
    let bwa_result = silence_stderr(|| aligner.align_fastq_records(&reads, true, true, threads))??;
    
    bwa_result
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


