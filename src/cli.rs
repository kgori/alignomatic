use anyhow::{anyhow, Result};
use clap::Parser;
#[allow(unused_imports)]
use log::{debug, info, warn};
use crate::utils::check_file_exists;

#[derive(Parser)]
pub struct CliOptions {
    #[arg(short = '1', long, value_name = "FILE")]
    pub fastq_first: std::path::PathBuf,

    #[arg(short = '2', long, value_name = "FILE")]
    pub fastq_second: std::path::PathBuf,

    #[arg(
        short = 'i',
        long,
        value_name = "FILES",
        value_delimiter = ',',
        help = "Reference files to map against. Must be Fasta format, and must have set of BWA index files."
    )]
    pub index: Vec<std::path::PathBuf>,

    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        help = "Output folder for all fastq files. Will be created if it doesn't exist.",
        required = true
    )]
    pub output_folder: std::path::PathBuf,

    #[arg(short, long, default_value = "66666")]
    pub batch_size: usize,

    #[arg(
        short,
        long,
        default_value = "1",
        help = "Number of threads to use. One thread will be reserved for the main program; any extra threads will be used for read mapping."
    )]
    pub threads: usize,

    #[arg(
        long,
        default_value = "30",
        help = "Minimum size of a block of bases that will be considered unmapped."
    )]
    pub min_block_size: usize,

    #[arg(
        long,
        default_value = "10",
        help = "Minimum average base quality of a block of bases that will be considered unmapped."
    )]
    pub min_block_quality: f32,
}

pub fn parse_cli() -> Result<CliOptions> {
    let mut opts = CliOptions::parse();

    check_file_exists(&opts.fastq_first)?;
    check_file_exists(&opts.fastq_second)?;
    opts.index
        .iter()
        .try_for_each(|reference_file| -> Result<()> {
            check_file_exists(reference_file)?;
            for bwa_extension in &["amb", "ann", "bwt", "pac", "sa"] {
                let index_file = std::path::PathBuf::from(format!(
                    "{}.{}",
                    reference_file.to_string_lossy(),
                    bwa_extension
                ));
                check_file_exists(&index_file)?;
            }
            Ok(())
        })?;

    if opts.batch_size < 1 {
        return Err(anyhow!("Batch size must be greater than 0"));
    }

    if opts.threads < 1 {
        return Err(anyhow!("Number of threads must be greater than 0"));
    }

    if opts.min_block_size < 1 {
        return Err(anyhow!("Minimum block size must be greater than 0"));
    }

    if opts.min_block_quality < 0.0 {
        return Err(anyhow!("Minimum block quality must be greater than 0"));
    }

    info!(target: "PARAMS", "First fastq file: {:?}", opts.fastq_first);
    info!(target: "PARAMS", "Second fastq file: {:?}", opts.fastq_second);
    info!(target: "PARAMS", "Index files:");
    for index in &opts.index {
        info!(target: "PARAMS", "{}", index.to_string_lossy());
    }
    info!(target: "PARAMS", "Batch size: {}", opts.batch_size);
    info!(target: "PARAMS", "Threads: {}", opts.threads);
    info!(target: "PARAMS", "Minimum block size: {}", opts.min_block_size);
    info!(target: "PARAMS", "Minimum block quality: {}", opts.min_block_quality);
    info!(target: "PARAMS", "Output folder: {:?}", opts.output_folder);

    if opts.threads > 1 {
        opts.threads -= 1;
    }

    Ok(opts)
}

