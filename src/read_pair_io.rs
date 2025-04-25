use crate::read_types::{MappedReadPair, ReadPair};
use crate::utils::bam_to_fastq;
use anyhow::{anyhow, Result};
use bgzf::Writer as BgzfWriter;
use bio::io::fastq::{self, FastqRead};
use log::debug;
use rust_htslib::bam;
use std::fs::File;
use std::io;

/// ReadPairIterator handles reading input from paired fastq files.
/// API summary:
/// new(file1, file2) - create a new ReadPairIterator from paths to two fastq
/// files.
/// next() - Iterator interface. Gets the next read pair.
/// batch(size) - returns a batch of reads with the first- and second-in-pair
/// interleaved.
/// The batch size is the number of read pairs to return (so the total number
/// of reads is double this).
pub struct ReadPairIterator {
    reader1: fastq::Reader<io::BufReader<Box<dyn io::Read>>>,
    reader2: fastq::Reader<io::BufReader<Box<dyn io::Read>>>,
}

impl ReadPairIterator {
    pub fn new(file1: std::path::PathBuf, file2: std::path::PathBuf) -> Result<Self> {
        debug!(target: "IO", "Opening files: {:?}, {:?}", file1, file2);
        let file1 = File::open(file1)?;
        let file2 = File::open(file2)?;

        let (reader1, format) = niffler::get_reader(Box::new(file1))?;
        let (reader2, format2) = niffler::get_reader(Box::new(file2))?;

        debug!(target: "IO", "Detected file1 compression: {:?}", format);
        debug!(target: "IO", "Detected file2 compression: {:?}", format2);

        let reader1 = fastq::Reader::new(reader1);
        let reader2 = fastq::Reader::new(reader2);
        Ok(Self { reader1, reader2 })
    }

    pub fn take_pairs(&mut self, batch_size: usize) -> Vec<ReadPair> {
        let mut batch = Vec::with_capacity(batch_size);

        for _ in 0..batch_size {
            if let Some(read_pair) = self.next() {
                batch.push(read_pair);
            } else {
                break;
            }
        }
        batch
    }

    pub fn take_bases(&mut self, base_pairs: usize) -> Vec<ReadPair> {
        let mut batch = Vec::new();
        let mut collected = 0;

        while collected < base_pairs {
            if let Some(read_pair) = self.next() {
                let bases = read_pair.read1.seq().len() + read_pair.read2.seq().len();
                batch.push(read_pair);
                collected += bases;
            } else {
                break;
            }
        }
        debug!(target: "IO", "Collected {} read pairs, {} bases", batch.len(), collected);
        batch
    }
}

impl Iterator for ReadPairIterator {
    type Item = ReadPair;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record1 = fastq::Record::new();
        let mut record2 = fastq::Record::new();

        let result1 = self.reader1.read(&mut record1);
        let result2 = self.reader2.read(&mut record2);

        if result1.is_ok() && result2.is_ok() {
            let read_pair = ReadPair::new(record1, record2);
            read_pair.ok()
        } else {
            None
        }
    }
}

pub type FastqWriter = fastq::Writer<BgzfWriter<std::fs::File>>;

pub fn create_bgzf_fastq_writer(path: &std::path::Path) -> Result<FastqWriter> {
    let file = std::fs::File::create(path)?;
    let bgzf_writer = BgzfWriter::new(file, 6.try_into()?);
    Ok(fastq::Writer::new(bgzf_writer))
}

/// Writes all alignments (primary, secondary, and supplementary) of a mapped read pair to the BAM file.
pub fn write_mapped_pair_to_bam(
    bam_writer: &mut bam::Writer,
    mapped_pair: &MappedReadPair,
) -> Result<usize> {
    let mut write_count = 0;

    for bam_record in mapped_pair.read1.iter() {
        bam_writer.write(bam_record)?;
        write_count += 1;
    }

    for bam_record in mapped_pair.read2.iter() {
        bam_writer.write(bam_record)?;
        write_count += 1;
    }

    Ok(write_count)
}

/// Writes only the primary alignment of a mapped read pair to the FASTQ files.
/// Errors if either read lacks a primary alignment.
pub fn write_mapped_pair_to_fastq(
    fastq_writer_1: &mut FastqWriter,
    fastq_writer_2: &mut FastqWriter,
    mapped_pair: &MappedReadPair,
) -> Result<()> {
    // Extract the primary alignment for read1
    let primary_read1: Option<&bam::Record> = mapped_pair
        .read1
        .iter()
        .find(|record| !record.is_secondary() && !record.is_supplementary());

    // Extract the primary alignment for read2
    let primary_read2: Option<&bam::Record> = mapped_pair
        .read2
        .iter()
        .find(|record| !record.is_secondary() && !record.is_supplementary());

    // Ensure both reads have primary alignments
    if primary_read1.is_none() || primary_read2.is_none() {
        return Err(anyhow!(
            "Primary alignment missing for one or both reads in pair: {}",
            mapped_pair.id()
        ));
    }

    // Convert primary alignments to FASTQ records
    let fastq_record1 = bam_to_fastq(primary_read1.unwrap())?;
    let fastq_record2 = bam_to_fastq(primary_read2.unwrap())?;

    // Write both FASTQ records
    fastq_writer_1.write_record(&fastq_record1)?;
    fastq_writer_2.write_record(&fastq_record2)?;

    Ok(())
}
