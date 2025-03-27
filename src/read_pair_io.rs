use crate::mapping_status::MappingStatus;
use anyhow::{anyhow, Result};
use bio::io::fastq::{self, FastqRead};
use log::debug;
use niffler;
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

    pub fn batch(&mut self, batch_size: usize) -> Vec<ReadPair> {
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
            if let Ok(read_pair) = read_pair {
                return Some(read_pair);
            } else {
                // error!("Failed to create a read pair: {:?}", read_pair.err());
                return None;
            }
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct ReadPair {
    pub id: String,
    pub read1: fastq::Record,
    pub read2: fastq::Record,
}

impl ReadPair {
    pub fn new(read1: fastq::Record, read2: fastq::Record) -> Result<Self> {
        let read1_id = Self::strip_comment(read1.id());
        if read1_id != Self::strip_comment(read2.id()) {
            return Err(anyhow!(
                "Read IDs do not match: {} != {}",
                read1.id(),
                read2.id()
            ));
        }
        if let Some(id) = read1_id {
            let id = id.to_string();
            let read_pair = Self { id, read1, read2 };
            return Ok(read_pair);
        }
        Err(anyhow!("Read IDs are empty"))
    }

    // private function to strip away comments from the read id
    fn strip_comment(id: &str) -> Option<&str> {
        id.split_whitespace().next()
    }
}

#[derive(Debug, Clone)]
pub struct SequencingRead {
    pub fastq: fastq::Record,
    pub alignments: Vec<bam::Record>,
    pub status: MappingStatus,
}

impl SequencingRead {
    pub fn new(fastq: fastq::Record) -> Self {
        SequencingRead {
            fastq,
            alignments: Vec::new(),
            status: MappingStatus::Unknown,
        }
    }
}

#[derive(Debug)]
pub struct MappedReadPair {
    #[allow(dead_code)]
    id: String,
    pub read1: Vec<bam::Record>,
    pub read2: Vec<bam::Record>,
}

impl MappedReadPair {
    pub fn new(id: &str) -> Self {
        MappedReadPair {
            id: id.to_string(),
            read1: Vec::new(),
            read2: Vec::new(),
        }
    }

    pub fn insert(&mut self, read: bam::Record) {
        if read.is_first_in_template() {
            self.read1.push(read);
        } else {
            self.read2.push(read);
        }
    }
}
