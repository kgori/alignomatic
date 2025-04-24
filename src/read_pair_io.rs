use crate::read_types::ReadPair;
use anyhow::Result;
use bio::io::fastq::{self, FastqRead};
use log::debug;
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
            if let Ok(read_pair) = read_pair {
                Some(read_pair)
            } else {
                // error!("Failed to create a read pair: {:?}", read_pair.err());
                None
            }
        } else {
            None
        }
    }
}
