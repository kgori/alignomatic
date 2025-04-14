use crate::utils::{check_file_exists, normalize_path};
use anyhow::{anyhow, Result};
use clap::Parser;
#[allow(unused_imports)]
use log::{debug, info, warn};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Parser)]
struct CliOptions {
    /// Optional path to a config file (JSON format)
    #[arg(long)]
    config: Option<PathBuf>,

    #[arg(short = '1', long, value_name = "FILE")]
    fastq_first: Option<PathBuf>,

    #[arg(short = '2', long, value_name = "FILE")]
    fastq_second: Option<PathBuf>,

    #[arg(
        short = 'i',
        long,
        value_name = "FILES",
        value_delimiter = ',',
        help = "Reference files to map against. Must be Fasta format, and must have set of BWA index files."
    )]
    index: Option<Vec<PathBuf>>,

    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        help = "Output folder for all fastq files. Will be created if it doesn't exist."
    )]
    output_folder: Option<PathBuf>,

    #[arg(short, long, help = "Batch size to process, in base pairs per thread. Default is 10000000.")]
    batch_size: Option<usize>,

    #[arg(
        short,
        long,
        help = "Number of threads to use. One thread will be reserved for the main program; any extra threads will be used for read mapping."
    )]
    threads: Option<usize>,

    #[arg(
        long,
        help = "Minimum size of a block of bases that will be considered unmapped."
    )]
    min_block_size: Option<usize>,

    #[arg(
        long,
        help = "Minimum average base quality of a block of bases that will be considered unmapped."
    )]
    min_block_quality: Option<f32>,
}

#[derive(Debug, Deserialize, Default)]
struct ConfigFileOptions {
    fastq_first: Option<PathBuf>,
    fastq_second: Option<PathBuf>,
    index: Option<Vec<PathBuf>>,
    output_folder: Option<PathBuf>,
    batch_size: Option<usize>,
    threads: Option<usize>,
    min_block_size: Option<usize>,
    min_block_quality: Option<f32>,
}

#[derive(Debug, Serialize)]
pub struct ProgramOptions {
    pub fastq_first: PathBuf,
    pub fastq_second: PathBuf,
    pub index: Vec<PathBuf>,
    pub output_folder: PathBuf,
    pub batch_size: usize,
    pub threads: usize,
    pub min_block_size: usize,
    pub min_block_quality: f32,
}

fn load_config(path: &PathBuf) -> ConfigFileOptions {
    let content = std::fs::read_to_string(path).expect("Failed to read config file");
    serde_json::from_str(&content).expect("Failed to parse config file")
}

fn merge_options(cli: CliOptions, config: ConfigFileOptions) -> Result<ProgramOptions> {
    let options = ProgramOptions {
        fastq_first: normalize_path(cli
            .fastq_first
            .or(config.fastq_first)
            .expect("No first fastq file provided"))?,
        fastq_second: normalize_path(cli
            .fastq_second
            .or(config.fastq_second)
            .expect("No second fastq file provided"))?,
        index: cli
            .index
            .or(config.index)
            .expect("No index files provided")
            .iter()
            .map(normalize_path)
            .collect::<Result<Vec<PathBuf>>>()?,
        output_folder: normalize_path(cli
            .output_folder
            .or(config.output_folder)
            .expect("No output folder provided"))?,
        batch_size: cli.batch_size.or(config.batch_size).unwrap_or(10_000_000),
        threads: cli.threads.or(config.threads).unwrap_or(1),
        min_block_size: cli.min_block_size.or(config.min_block_size).unwrap_or(30),
        min_block_quality: cli
            .min_block_quality
            .or(config.min_block_quality)
            .unwrap_or(10.0),
    };
    Ok(options)
}

fn write_config(config: &ProgramOptions) -> Result<()> {
    let json_str = serde_json::to_string_pretty(config)?;
    let path = config.output_folder.join("config").with_extension("json");
    std::fs::write(&path, json_str)?;
    info!(target: "PARAMS", "Wrote config file to {:?}", path);
    Ok(())
}

fn check_options(opts: &ProgramOptions) -> Result<()> {
    check_file_exists(&opts.fastq_first)?;
    check_file_exists(&opts.fastq_second)?;
    opts.index
        .iter()
        .try_for_each(|reference_file| -> Result<()> {
            check_file_exists(reference_file)?;
            for bwa_extension in &["amb", "ann", "bwt", "pac", "sa"] {
                let index_file = PathBuf::from(format!(
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

    Ok(())
}

fn create_output_folder(opts: &ProgramOptions) -> Result<()> {
    if opts.output_folder.exists() {
        warn!(target: "IO", "Output folder {} already exists. Files may be overwritten.", opts.output_folder.display());
    } else {
        info!(target: "IO", "Creating output folder {}", opts.output_folder.display());
        std::fs::create_dir(&opts.output_folder)?;
    }

    let workspace = opts.output_folder.join("workspace");
    let results = opts.output_folder.join("results");
    std::fs::create_dir_all(workspace)?;
    std::fs::create_dir_all(results)?;
    Ok(())
}

pub fn get_program_options() -> Result<ProgramOptions> {
    let cli_opts = CliOptions::parse();
    let config_opts = cli_opts
        .config
        .as_ref()
        .map(load_config)
        .unwrap_or_default();
    let final_opts = merge_options(cli_opts, config_opts)?;
    check_options(&final_opts)?;
    create_output_folder(&final_opts)?;
    write_config(&final_opts)?;
    Ok(final_opts)
}
