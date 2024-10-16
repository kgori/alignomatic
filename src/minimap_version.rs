#![allow(unused_variables)]
#![allow(unused_imports)]
#![allow(dead_code)]
use clap::Parser;
use bio::io::fastq::{self, FastqRead};
use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::BufReader;
use flate2::read::GzDecoder;

use rust_htslib::bam;
use rust_htslib::bam::{Format, Header, HeaderView};


#[derive(Parser)]
struct Cli {
    #[arg(short = '1', long, value_name = "FILE")]
    fastq_first: std::path::PathBuf,

    #[arg(short = '2', long, value_name = "FILE")]
    fastq_second: std::path::PathBuf,
    
    #[arg(short = 'i', long, value_name = "FILE")]
    index: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opts = Cli::parse();
    println!("First fastq file: {:?}", opts.fastq_first);
    println!("Second fastq file: {:?}", opts.fastq_second);
    println!("Index file: {:?}", opts.index);
    do_work(&opts)?;
    Ok(())
}

fn create_minimap_aligner_from_index(index: &std::path::PathBuf) -> Result<minimap2::Aligner> {
    let aligner = minimap2::Aligner::builder()
        .with_cigar()
        .with_index(index, None)
        .map_err(|err| anyhow!("Failed to create minimap2 aligner: {}", err))?;
    Ok(aligner)
}

pub struct ReadPairIterator {
    reader1: fastq::Reader<BufReader<GzDecoder<File>>>,
    reader2: fastq::Reader<BufReader<GzDecoder<File>>>,
}

impl ReadPairIterator {
    pub fn new(file1: std::path::PathBuf, file2: std::path::PathBuf) -> Result<Self> {
        let file1 = File::open(file1)?;
        let file2 = File::open(file2)?;
        let reader1 = fastq::Reader::new(GzDecoder::new(file1));
        let reader2 = fastq::Reader::new(GzDecoder::new(file2));
        Ok(Self { reader1, reader2 })
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
            Some(ReadPair::new(record1, record2).unwrap())
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct ReadPair {
    id: String,
    read1: fastq::Record,
    read2: fastq::Record,
}

impl ReadPair {
    pub fn new(read1: fastq::Record, read2: fastq::Record) -> Result<Self> {
        if Self::strip_comment(read1.id()) != Self::strip_comment(read2.id()) {
            return Err(anyhow!("Read IDs do not match: {} != {}", read1.id(), read2.id()));
        }
        let id = Self::strip_comment(read1.id()).to_string();
        Ok(Self { id, read1, read2 })
    }

    // private function to strip away comments from the read id
    fn strip_comment(id: &str) -> &str {
        id.split_whitespace().next().unwrap()
    }
}

fn do_work(opts: &Cli) -> Result<()> {
    let aligner = create_minimap_aligner_from_index(&opts.index)?;
    let mut header = Header::new();
    aligner.populate_header(&mut header);
    let header_view = HeaderView::from_header(&header);


    let mut reader = ReadPairIterator::new(opts.fastq_first.clone(), opts.fastq_second.clone())?;
    let read_pair = reader.next().unwrap();
    println!("{:?}", read_pair);
    let read_pair = reader.next().unwrap();
    println!("{:?}", read_pair);
    let hits = aligner.map_to_sam(read_pair.read1.seq(), Some(read_pair.read1.qual()), Some(read_pair.id.as_bytes()), &header_view, None, None);
    let mut sam_writer = bam::Writer::from_stdout(&header, Format::Sam).unwrap();
    if let Ok(hits) = hits {
        for record in hits {
            sam_writer.write(&record).unwrap();
        }
    }
    Ok(())
}
