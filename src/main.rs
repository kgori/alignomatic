#![allow(dead_code)]
#![allow(unused_imports)]

use anyhow::Result;
use bio::io::fastq;
use log::{debug, error, info, warn};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::BTreeMap;
use std::path::PathBuf;

mod aligner;
mod cli;
mod mapping_status;
mod output;
mod read_pair_io;
mod utils;

use aligner::Aligner;
use mapping_status::{get_mapping_status, MappingStatus, MappingStatus::*};
use read_pair_io::{MappedReadPair, ReadPairIterator};
use utils::{bam_to_fastq, create_bgzf_fastq_writer};

fn main() -> Result<()> {
    env_logger::init();

    let opts = cli::parse_cli()?;
    if opts.output_folder.exists() {
        warn!(target: "IO", "Output folder {} already exists. Files may be overwritten.", opts.output_folder.display());
    } else {
        info!(target: "IO", "Creating output folder {}", opts.output_folder.display());
        std::fs::create_dir(&opts.output_folder)?;
    }
    std::fs::create_dir_all(&opts.output_folder.join("workspace"))?;
    std::fs::create_dir_all(&opts.output_folder.join("results"))?;
    let bam_files = generate_alignments(&opts)?;
    let result = post_process_alignments(&bam_files, &opts);
    if result.is_ok() {
        info!(target: "IO", "Finished processing alignments");
    } else {
        error!(target: "IO", "Error processing alignments");
    }
    Ok(())
}

fn generate_alignments(opts: &cli::CliOptions) -> Result<Vec<PathBuf>> {
    info!(target: "IO", "Accessing Fastq read pairs from disk");

    let mut bam_files = Vec::new();
    let mut fastq1 = opts.fastq_first.clone();
    let mut fastq2 = opts.fastq_second.clone();

    for index in &opts.index {
        info!(target: "BWA", "Loading BWA index from {}", index.display());
        let mut read_pair_iter = ReadPairIterator::new(fastq1.clone(), fastq2.clone())?;

        let mut aligner = Aligner::new(&index, &opts.output_folder.join("workspace"))?;

        loop {
            let read_pairs = read_pair_iter.batch(opts.batch_size);
            let n_read_pairs = read_pairs.len();

            if n_read_pairs == 0 {
                break;
            } else {
                info!(target: "BWA", "Aligning {} read pairs", read_pairs.len());
                aligner.process_batch(&read_pairs, &opts)?;
            }
        }
        bam_files.push(aligner.bamfile().clone());
        (fastq1, fastq2) = aligner.fastqfiles();
    }

    Ok(bam_files)
}

fn post_process_alignments(bam_files: &[PathBuf], opts: &cli::CliOptions) -> Result<()> {
    let mut alignments: std::collections::BTreeMap<String, MappedReadPair> =
        std::collections::BTreeMap::new();

    for bam_file in bam_files {
        let mut bam = bam::Reader::from_path(bam_file)?;
        for record in bam.records() {
            let record = record?;
            let read_name = String::from_utf8(record.qname().to_vec())?;
            let mapped_pair = alignments
                .entry(read_name.clone())
                .or_insert(MappedReadPair::new(&read_name));
            mapped_pair.insert(record);
        }
    }

    let mut unmapped_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_uu.1.fq.gz"))?;
    let mut unmapped_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_uu.2.fq.gz"))?;
    let mut mapped_unmapped_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_mu.1.fq.gz"))?;
    let mut mapped_unmapped_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_mu.2.fq.gz"))?;
    let mut unmapped_mapped_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_um.1.fq.gz"))?;
    let mut unmapped_mapped_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_um.2.fq.gz"))?;
    let mut fragmentary_unmapped_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_fu.1.fq.gz"))?;
    let mut fragmentary_unmapped_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_fu.2.fq.gz"))?;
    let mut unmapped_fragmentary_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_uf.1.fq.gz"))?;
    let mut unmapped_fragmentary_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_uf.2.fq.gz"))?;
    let mut fragmentary_mapped_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_fm.1.fq.gz"))?;
    let mut fragmentary_mapped_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_fm.2.fq.gz"))?;
    let mut mapped_fragmentary_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_mf.1.fq.gz"))?;
    let mut mapped_fragmentary_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_mf.2.fq.gz"))?;

    let mut fragmentary_pairs1 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_ff.1.fq.gz"))?;
    let mut fragmentary_pairs2 =
        create_bgzf_fastq_writer(&opts.output_folder.join("results/reads_ff.2.fq.gz"))?;

    for (id, read_pair) in alignments {
        let status1 = get_mapping_status(&read_pair.read1, &opts)?;
        let status2 = get_mapping_status(&read_pair.read2, &opts)?;
        match (status1, status2) {
            (Mapped, Mapped) => {
                ();
            }
            (Unmapped, Unmapped) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                unmapped_pairs1.write_record(&f1)?;
                unmapped_pairs2.write_record(&f2)?;
            }
            (Unmapped, Mapped) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                unmapped_mapped_pairs1.write_record(&f1)?;
                unmapped_mapped_pairs2.write_record(&f2)?;
            }
            (Mapped, Unmapped) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                mapped_unmapped_pairs1.write_record(&f1)?;
                mapped_unmapped_pairs2.write_record(&f2)?;
            }
            (Fragmentary(_), Unmapped) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                fragmentary_unmapped_pairs1.write_record(&f1)?;
                fragmentary_unmapped_pairs2.write_record(&f2)?;
            }
            (Unmapped, Fragmentary(_)) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                unmapped_fragmentary_pairs1.write_record(&f1)?;
                unmapped_fragmentary_pairs2.write_record(&f2)?;
            }
            (Fragmentary(_), Mapped) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                fragmentary_mapped_pairs1.write_record(&f1)?;
                fragmentary_mapped_pairs2.write_record(&f2)?;
            }
            (Mapped, Fragmentary(_)) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                mapped_fragmentary_pairs1.write_record(&f1)?;
                mapped_fragmentary_pairs2.write_record(&f2)?;
            }
            (Fragmentary(_), Fragmentary(_)) => {
                let f1 = convert_bam_to_fastq(&read_pair.read1)?;
                let f2 = convert_bam_to_fastq(&read_pair.read2)?;
                fragmentary_pairs1.write_record(&f1)?;
                fragmentary_pairs2.write_record(&f2)?;
            }
            (Suspicious, _) | (_, Suspicious) => {
                error!("Suspicious read mapping: {}", id);
            }
            _ => {
                error!("Unknown read mapping: {}", id);
            }
        }
    }

    Ok(())
}

fn convert_bam_to_fastq(alignments: &[bam::Record]) -> Result<fastq::Record> {
    let mut fastqs = Vec::new();
    for alignment in alignments.iter() {
        if alignment.is_secondary() || alignment.is_supplementary() {
            continue;
        }
        fastqs.push(bam_to_fastq(alignment.clone())?);
    }
    match fastqs.len() {
        0 => Err(anyhow::anyhow!("No primary alignments")),
        1 => Ok(fastqs[0].clone()),
        2.. => {
            let read1_id = fastqs[0].id();
            let read1_seq = fastqs[0].seq();
            let read1_qual = fastqs[0].qual();
            for fastq in fastqs.iter().skip(1) {
                if read1_id != fastq.id() {
                    return Err(anyhow::anyhow!("Read IDs do not match"));
                }
                if read1_seq != fastq.seq() {
                    return Err(anyhow::anyhow!("Read sequences do not match"));
                }
                if read1_qual != fastq.qual() {
                    return Err(anyhow::anyhow!("Read qualities do not match"));
                }
            }
            Ok(fastqs[0].clone())
        }
    }
}
