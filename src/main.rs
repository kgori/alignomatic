#![allow(dead_code)]

use anyhow::Result;
use bam_io::CollatedBamReader;
use bio::io::fastq;
use log::{debug, error, info, warn};
use rust_htslib::bam::{self, Read};
use std::path::PathBuf;

mod aligner;
mod bam_io;
mod checkpoint;
mod cli;
mod mapping_status;
mod read_pair_io;
mod read_types;
mod utils;

use aligner::Aligner;
use mapping_status::{get_mapping_status, MappingStatus::*};
use read_types::MappedReadPair;
use read_pair_io::{create_bgzf_fastq_writer, write_mapped_pair_to_bam, ReadPairIterator};
use utils::{bam_to_fastq, extract_primary_read, read_pair_is_unmapped};

fn main() -> Result<()> {
    if std::env::args().len() == 1 {
        println!("No arguments provided. Use --help for usage information.");
        std::process::exit(0);
    }

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let opts: cli::ProgramOptions = cli::get_program_options()?;

    let bam_files = if opts.bam_input.is_some() {
        generate_alignments_from_bam(&opts)?
    } else {
        generate_alignments(&opts)?
    };
    let result = post_process_alignments(&bam_files, &opts, 100_000);
    if result.is_ok() {
        info!(target: "IO", "Finished processing alignments");
    } else {
        error!(target: "IO", "Error processing alignments");
    }

    Ok(())
}

fn get_bam_header_from_bam(bamfile: &PathBuf) -> Result<bam::Header> {
    let bam_reader = bam::Reader::from_path(bamfile)?;
    let header = bam::Header::from_template(bam_reader.header());
    Ok(header)
}

fn generate_alignments_from_bam(opts: &cli::ProgramOptions) -> Result<Vec<PathBuf>> {
    info!(target: "IO", "Accessing input BAM file from disk");

    let workspace = opts.output_folder.join("workspace");
    let mut bam_files = Vec::new();

    let bam_header = get_bam_header_from_bam(&opts.bam_input.clone().unwrap())?;
    let bam_writer_filename = workspace.join("initial_output.bam");
    let mut bam_writer = bam::Writer::from_path(
        bam_writer_filename.clone(),
        &bam_header,
        bam::Format::Bam,
    )?;

    let fastqout1_filename = workspace.join("initial_output.1.fq.gz");
    let mut fastqout1 =
        create_bgzf_fastq_writer(&fastqout1_filename)?;

    let fastqout2_filename = workspace.join("initial_output.2.fq.gz");
    let mut fastqout2 =
        create_bgzf_fastq_writer(&fastqout2_filename)?;

    info!(target: "IO", "Collating BAM file {}", opts.bam_input.clone().unwrap().display());
    let bam_reader = CollatedBamReader::new(
        opts.bam_input.clone().unwrap().as_path(),
        Some(workspace.as_path()),
        Some(opts.threads),
    )?;

    let mut bam_write_count = 0;
    let mut fastq_write_count: usize = 0;
    let mut seen_count: usize = 0;
    
    for mapped_pair in bam_reader {
        seen_count += 1;
        let mapped_pair = mapped_pair?;
        if mapped_pair.read1.is_empty() || mapped_pair.read2.is_empty() {
            warn!(target: "IO",
                "Read pair {} is missing at least one of its pair. Skipping.",
                mapped_pair.id());
            continue;
        }
        if read_pair_is_unmapped(&mapped_pair, opts)? {
            bam_write_count += write_mapped_pair_to_bam(&mut bam_writer, &mapped_pair)?;

            // Find the primary alignment in mapped_pair.read1 and convert it to fastq
            let primary_read1 = extract_primary_read(&mapped_pair.read1)?;
            fastqout1.write_record(&primary_read1)?;

            let primary_read2 = extract_primary_read(&mapped_pair.read2)?;
            fastqout2.write_record(&primary_read2)?;

            fastq_write_count += 1;

            if fastq_write_count % 1_000_000 == 0 {
                debug!(target: "IO",
                    "Wrote {} reads to fastq files",
                    fastq_write_count);
            }
        }
        if seen_count % 1_000_000 == 0 {
            debug!(target: "IO",
                "Processed {} read pairs",
                seen_count);
        }
    }

    debug!(target: "IO",
        "Wrote {} reads to BAM file and {} reads to fastq files",
        bam_write_count,
        fastq_write_count);

    drop(bam_writer);

    bam_files.push(bam_writer_filename);

    process_with_index(
        opts,
        &mut bam_files,
        fastqout1_filename,
        fastqout2_filename,
         workspace)?;
    Ok(bam_files)
}

fn generate_alignments(opts: &cli::ProgramOptions) -> Result<Vec<PathBuf>> {
    info!(target: "IO", "Accessing Fastq read pairs from disk");
    let fastq1 = opts.fastq_first.clone().unwrap();
    let fastq2 = opts.fastq_second.clone().unwrap();

    let mut bam_files = Vec::new();

    let workspace = opts.output_folder.join("workspace");

    process_with_index(opts, &mut bam_files, fastq1, fastq2, workspace)?;

    Ok(bam_files)
}

fn process_with_index(
    opts: &cli::ProgramOptions,
    bam_files: &mut Vec<PathBuf>,
    mut fastq1: PathBuf,
    mut fastq2: PathBuf,
    workspace: PathBuf
) -> Result<()> {
    for index in &opts.index {
        info!(target: "Aligner", "Loading index {}", index.display());

        // CHECKPOINTING:
        // Look for presence of a checkpoint file
        //   if found:
        //     1 - compute the checkpoint hash for the current inputs (ref infq1 infq2 bam outfq1 outfq2)
        //     2 - compare the checkpoint hash with the one in the checkpoint file
        //     3 - if they match, skip the loop and continue with the next index
        //   otherwise:
        //     1 - Do checkpoint hashing work in background
        //     2 - run the loop
        //     3 - write the checkpoint file with the hash data

        info!(target: "Checkpoint", "Looking for checkpoint file for index {}", index.display());
        let (bamfile, fqout1, fqout2, checkpoint_file) =
            aligner::output_filenames(index, &workspace)?;
        debug!(target: "Checkpoint", "Checkpoint file: {}", checkpoint_file.display());
        debug!(target: "Checkpoint", "Output bam file: {}", bamfile.display());
        debug!(target: "Checkpoint", "Output fastq file 1: {}", fqout1.display());
        debug!(target: "Checkpoint", "Output fastq file 2: {}", fqout2.display());

        if checkpoint_file.exists() {
            info!(target: "Checkpoint",
                "Checkpoint file found for index {}",
                index.display());

            let checkpoint_on_disk = checkpoint::read_checkpoint(&checkpoint_file)?;
            let checkpoint_computed = checkpoint::Checkpoint::create(
                index, &fastq1, &fastq2, &bamfile, &fqout1, &fqout2, opts.batch_size,
                opts.min_block_size, opts.min_block_quality,
            )?;

            if checkpoint_on_disk.matches(&checkpoint_computed) {
                info!(target: "Checkpoint", "Checkpoint verified for index {}. Skipping alignment.", index.display());
                bam_files.push(bamfile);
                (fastq1, fastq2) = (fqout1, fqout2);
                continue;
            } else {
                info!(target: "Checkpoint",
                    "Checkpoint mismatch for index {}. Recomputing alignment.",
                    index.display());
                debug!(target: "Checkpoint",
                    "Checkpoint on disk: {:?}",
                    checkpoint_on_disk);
                debug!(target: "Checkpoint",
                    "Checkpoint computed: {:?}",
                    checkpoint_computed);
            }
        } else {
            info!(target: "Checkpoint",
                "Checkpoint not found for index {}. Computing alignment.",
                index.display());
        }

        // If we got here it means there is work to be done (no checkpoint file or hash mismatch)
        // Create the aligner
        {
            let mut aligner = Aligner::new(
                index, &workspace,
                Some(opts.threads), // Writing threads: use multiple threads for writing and aligning
            )?;
            info!(target: "Aligner", "Aligning reads to index {}", index.display());
            let mut read_pair_iter = ReadPairIterator::new(fastq1.clone(), fastq2.clone())?;
            assert_eq!(aligner.bamfile(), &bamfile);
            assert_eq!(aligner.fastqfiles(), (&fqout1, &fqout2));
            let mut n_readpairs_processed: usize = 0;

            // Process all the reads
            loop {
                let read_pairs = read_pair_iter.take_bases(opts.threads * opts.batch_size);
                let n_read_pairs = read_pairs.len();

                if n_read_pairs == 0 {
                    break;
                } else {
                    info!(target: "BWA",
                        "Aligning {} read pairs",
                        read_pairs.len());
                    aligner.process_batch(&read_pairs, opts)?;
                    n_readpairs_processed += n_read_pairs;
                    info!(target: "BWA",
                        "{} read pairs processed",
                        n_readpairs_processed);
                }
            }
        }

        // Write checkpoint file
        info!(target: "Checkpoint",
            "Successfully completed. Writing checkpoint file for index {}",
            index.display());
        let ckpt =
            checkpoint::Checkpoint::create(index, &fastq1, &fastq2, &bamfile, &fqout1, &fqout2,
                opts.batch_size, opts.min_block_size, opts.min_block_quality)?;
        checkpoint::write_checkpoint(&ckpt, &checkpoint_file)?;
        bam_files.push(bamfile);
        (fastq1, fastq2) = (fqout1, fqout2);
    }

    Ok(())
}

fn post_process_alignments(bam_files: &[PathBuf], opts: &cli::ProgramOptions, batch_size: usize) -> Result<()> {
    let final_bamfile = bam_files.last().unwrap();
    let mut final_bam =
        bam_io::BufferedBamReader::new(bam::Reader::from_path(final_bamfile).unwrap());

    let mut bam_readers = bam_files
        .iter()
        .take(bam_files.len() - 1)
        .map(|bamfile| bam_io::BufferedBamReader::new(bam::Reader::from_path(bamfile).unwrap()))
        .collect::<Vec<_>>();

    info!(target: "Results", "Writing results to fastq files in {}", opts.output_folder.join("results").display());
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

    debug!(target: "Results", "Reading BAMs in batches of size {}", batch_size);
    while let Ok(batch) = final_bam.take_n_qnames(batch_size) {
        if batch.is_empty() {
            break;
        }

        let mut alignments: std::collections::BTreeMap<String, MappedReadPair> =
            std::collections::BTreeMap::new();
        let end_of_batch = String::from_utf8(batch.last().unwrap().qname().to_vec())?;

        for record in batch {
            let read_name = String::from_utf8(record.qname().to_vec())?;
            let mapped_pair = alignments
                .entry(read_name.clone())
                .or_insert(MappedReadPair::new(&read_name));
            mapped_pair.insert(record)?;
        }

        for bam_reader in bam_readers.iter_mut() {
            // This is a run-time check. All read pairs from the final bam file should be
            // present in all other bam files. If not, it's an error. This checks both read1
            // and read2.
            let mut check_map_read1: std::collections::BTreeMap<String, bool> =
                std::collections::BTreeMap::new();
            let mut check_map_read2: std::collections::BTreeMap<String, bool> =
                std::collections::BTreeMap::new();
            for key in alignments.keys() {
                check_map_read1.insert(key.clone(), false);
                check_map_read2.insert(key.clone(), false);
            }

            let batch = bam_reader.take_until_qname(&end_of_batch)?;
            for record in batch {
                let read_name = String::from_utf8(record.qname().to_vec())?;
                if let Some(mapped_pair) = alignments.get_mut(&read_name) {
                    if record.is_first_in_template() {
                        check_map_read1.insert(read_name, true); // Mark as found
                    } else {
                        check_map_read2.insert(read_name, true); // Mark as found
                    }
                    mapped_pair.insert(record)?;
                }
            }

            // Every key should be found, no exceptions.
            for (key, value) in check_map_read1.iter() {
                if !value {
                    return Err(anyhow::anyhow!(
                        "Read 1 of pair {} not found in all bam files",
                        key
                    ));
                }
            }

            for (key, value) in check_map_read2.iter() {
                if !value {
                    return Err(anyhow::anyhow!(
                        "Read 2 of pair {} not found in all bam files",
                        key
                    ));
                }
            }
        }

        debug!(target: "Results", "Processing {} read pairs", alignments.len());
        for (id, read_pair) in alignments {
            let status1 = get_mapping_status(&read_pair.read1, opts)?;
            let status2 = get_mapping_status(&read_pair.read2, opts)?;
            match (status1, status2) {
                (Mapped, Mapped) => { }
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
    }

    Ok(())
}

fn convert_bam_to_fastq(alignments: &[bam::Record]) -> Result<fastq::Record> {
    let mut fastqs = Vec::new();
    for alignment in alignments.iter() {
        if alignment.is_secondary() || alignment.is_supplementary() {
            continue;
        }
        fastqs.push(bam_to_fastq(alignment)?);
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
