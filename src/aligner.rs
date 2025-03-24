use crate::cli;
use crate::mapping_status::get_mapping_status;
use crate::mapping_status::MappingStatus::*;
use crate::read_pair_io::{MappedReadPair, ReadPair};
use crate::utils::{
    check_directory_exists, check_file_exists, create_bgzf_fastq_writer, silence_stderr,
    FastqWriter,
};
use anyhow::{anyhow, Result};
use bio::io::fastq;
use bwa::BwaAligner;
use log::{debug, warn};
use rust_htslib::bam;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A struct that represents a BWA aligner and its associated
/// files. It keeps track of the reference fasta file used for
/// the index, and a BAM file for writing output.
pub struct Aligner {
    aligner: BwaAligner,
    reference: PathBuf,
    output_dir: PathBuf,
    bam_filename: PathBuf,
    fastq_filename_1: PathBuf,
    fastq_filename_2: PathBuf,
    bam_writer: bam::Writer,
    fastq_writer_1: FastqWriter,
    fastq_writer_2: FastqWriter,
}

fn process_filename(filename: &PathBuf) -> Result<PathBuf> {
    let binding = filename.canonicalize()?;
    let processed = binding
        .file_name()
        .and_then(|s| s.to_str())
        .and_then(|s| Path::new(s).file_stem());

    match processed {
        Some(f) => Ok(PathBuf::from(f)),
        None => Err(anyhow!("Error processing filename {}", filename.display())),
    }
}

impl Aligner {
    pub fn new(reference: &PathBuf, output_dir: &PathBuf, writing_threads: Option<usize>) -> Result<Self> {
        check_file_exists(&reference)?;
        debug!(target: "Aligner", "Reference file exists: {}", reference.display());
        check_directory_exists(&output_dir)?;
        debug!(target: "Aligner", "Output directory exists: {}", output_dir.display());
        let aligner = silence_stderr(|| BwaAligner::from_path(&reference))??;
        let header = aligner.create_bam_header();
        let output_filename = process_filename(&reference)?;
        let bam_filename = output_dir
            .join(output_filename.clone())
            .with_extension("bam");
        let mut bam_writer = bam::Writer::from_path(&bam_filename, &header, bam::Format::Bam)?;
        if let Some(threads) = writing_threads {
            bam_writer.set_threads(threads)?;
        }
        let fastq_filename_1 = output_dir
            .join(output_filename.clone())
            .with_extension("1.fastq.gz");
        let fastq_filename_2 = output_dir
            .join(output_filename)
            .with_extension("2.fastq.gz");
        let fastq_writer_1 = create_bgzf_fastq_writer(&fastq_filename_1)?;
        let fastq_writer_2 = create_bgzf_fastq_writer(&fastq_filename_2)?;

        Ok(Self {
            aligner,
            reference: reference.clone(),
            output_dir: output_dir.clone(),
            bam_filename,
            bam_writer,
            fastq_filename_1,
            fastq_filename_2,
            fastq_writer_1,
            fastq_writer_2,
        })
    }

    pub fn bamfile(&self) -> PathBuf {
        self.bam_filename.clone()
    }

    pub fn fastqfiles(&self) -> (PathBuf, PathBuf) {
        (self.fastq_filename_1.clone(), self.fastq_filename_2.clone())
    }
    
    pub fn process_batch(&mut self, read_pairs: &[ReadPair], opts: &cli::ProgramOptions) -> Result<()> {
        let mut reads = Vec::<fastq::Record>::new();

        for read_pair in read_pairs.iter() {
            reads.push(read_pair.read1.clone());
            reads.push(read_pair.read2.clone());
        }

        let alignments = self.align_reads(&reads, opts.threads)?;

        let mut bam_write_count: usize = 0;
        let mut fastq_write_count: usize = 0;
        for input_pair in read_pairs {
            let id = &input_pair.id;
            if alignments.contains_key(id) {
                let mapped_pair = alignments.get(id).unwrap();
                let read1_status = get_mapping_status(&mapped_pair.read1, &opts)?;
                let read2_status = get_mapping_status(&mapped_pair.read2, &opts)?;
                match (read1_status, read2_status) {
                    (Mapped, Mapped) => {
                        continue;
                    }
                    (Suspicious, _) | (_, Suspicious) => {
                        warn!("Read pair {} has a suspicious alignment", input_pair.id);
                        continue;
                    }
                    _ => {
                        for bam_record in mapped_pair.read1.iter() {
                            self.bam_writer.write(&bam_record)?;
                            bam_write_count += 1;
                        }
                        for bam_record in mapped_pair.read2.iter() {
                            self.bam_writer.write(&bam_record)?;
                            bam_write_count += 1;
                        }
                        self.fastq_writer_1.write_record(&input_pair.read1)?;
                        self.fastq_writer_2.write_record(&input_pair.read2)?;
                        fastq_write_count += 1;
                    }
                }
            }
        }
        debug!(target: "IO", "Wrote {} alignments to {}", bam_write_count, self.bam_filename.display());
        debug!(target: "IO", "Wrote {} reads to {}", fastq_write_count, self.fastq_filename_1.display());
        debug!(target: "IO", "Wrote {} reads to {}", fastq_write_count, self.fastq_filename_2.display());

        Ok(())
    }

    fn align_reads(
        &self,
        reads: &[fastq::Record],
        threads: usize,
    ) -> Result<BTreeMap<String, MappedReadPair>> {
        debug!(target: "BWA", "Aligning {} reads", reads.len());
        let mut alignments = BTreeMap::new();
        let bwa_result = silence_stderr(|| {
            self.aligner
                .align_fastq_records(&reads, true, true, threads)
        })??;

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

    fn get_reference_filename(&self) -> Result<String> {
        // get the reference filename without its full path (just the basename)
        let reference_filename = self.reference.file_name().and_then(|s| s.to_str());
        match reference_filename {
            Some(f) => Ok(f.to_string()),
            None => Err(anyhow!("Error processing reference filename")),
        }
    }
}
