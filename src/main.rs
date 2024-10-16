#![allow(unused_variables)]
#![allow(dead_code)]
use anyhow::{anyhow, Result};
use clap::Parser;

use rust_htslib::bam;
use rust_htslib::bam::Format;

use bwa::BwaAligner;
use log::{info, warn};

mod read_pair_io;
use read_pair_io::ReadPairIterator;

mod utils;
use utils::silence_stderr;

#[derive(Parser)]
struct Cli {
    #[arg(short = '1', long, value_name = "FILE")]
    fastq_first: std::path::PathBuf,

    #[arg(short = '2', long, value_name = "FILE")]
    fastq_second: std::path::PathBuf,

    #[arg(
        short = 'i',
        long,
        value_name = "FILES",
        value_delimiter = ',',
        help = "Reference files to map against. Must be Fasta format, and must have set of BWA index files."
    )]
    index: Vec<std::path::PathBuf>,

    #[arg(short, long, default_value = "5000")]
    batch_size: usize,
}

fn main() -> Result<()> {
    env_logger::init();

    let opts = parse_cli()?;
    do_work(&opts)?;
    Ok(())
}

fn parse_cli() -> Result<Cli> {
    let opts = Cli::parse();

    check_file_exists(&opts.fastq_first)?;
    check_file_exists(&opts.fastq_second)?;
    opts.index.iter().try_for_each(|index| check_file_exists(index))?;
    
    info!(target: "Command line", "First fastq file: {:?}", opts.fastq_first);
    info!(target: "Command line", "Second fastq file: {:?}", opts.fastq_second);
    info!(target: "Command line", "Index files: {:?}", opts.index);
    info!(target: "Command line", "Batch size: {}", opts.batch_size);
    Ok(opts)
}

fn check_file_exists(file: &std::path::PathBuf) -> Result<()> {
    let existence = file.try_exists();
    match existence {
        Ok(true) => Ok(()),
        _ => Err(anyhow!("File not found: {}", file.to_string_lossy())),
    }
}

fn load_aligners(paths: &[std::path::PathBuf]) -> Result<Vec<BwaAligner>> {
    let mut aligners = Vec::<BwaAligner>::new();
    for pathbuf in paths.iter() {
        let aligner = silence_stderr(|| BwaAligner::from_path(pathbuf))??;
        aligners.push(aligner);
    };
    Ok(aligners)
}

fn do_work(opts: &Cli) -> Result<()> {
    info!(target: "IO", "Loading BWA indices and creating aligners");
    let aligners = load_aligners(&opts.index)?;
    let aligner = &aligners[0];
    let header = aligner.create_bam_header();

    info!(target: "IO", "Accessing Fastq read pairs from disk");
    let mut read_pair_iter =
        ReadPairIterator::new(opts.fastq_first.clone(), opts.fastq_second.clone())?;

    let mut writer = bam::Writer::from_stdout(&header, Format::Sam)?;
    let reads = read_pair_iter.batch(opts.batch_size);
    info!(target: "Aligner", "Aligning {} read pairs", reads.len());
    if reads.len() > 0 {
        let alignments = aligner.align_fastq_records(reads.as_slice(), true, true, 1)?;
        for aln in alignments.iter() {
            writer.write(aln)?;
        }
    } else {
        warn!(target: "IO", "No read pair found");
    }
    Ok(())
}
