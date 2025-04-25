use crate::mapping_status::MappingStatus;
use anyhow::{anyhow, Result};
use bio::io::fastq;
use rust_htslib::bam;

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

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn insert(&mut self, read: bam::Record) {
        if read.is_first_in_template() {
            self.read1.push(read);
        } else {
            self.read2.push(read);
        }
    }
}
