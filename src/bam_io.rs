use anyhow::{anyhow, Result};
use rust_htslib::bam::{self, Read};

pub struct BufferedBamReader {
    reader: bam::Reader,
    read: Option<bam::Record>,
}

impl BufferedBamReader {
    pub fn new(reader: bam::Reader) -> Self {
        Self { reader, read: None }
    }

    /// Reads records until `n` distinct qnames have been seen. Reads until the end of the nth group
    /// of records with the same qname. It's expected that the reads will be collated by qname, as will
    /// be the case for the bams written by this program to the workspace directory.
    pub fn take_n_qnames(&mut self, n: usize) -> Result<Vec<bam::Record>> {
        let mut batch = Vec::new();
        let mut record = bam::Record::new();
        let mut names_seen = 0;
        let mut prev_name: Option<String> = None;

        // Include any cached read if it exists
        if let Some(cached_record) = self.read.take() {
            batch.push(cached_record.clone());
            prev_name = Some(String::from_utf8(cached_record.qname().to_vec())?);
            names_seen += 1;
        }

        // Process remaining records
        while names_seen <= n {
            let result = self.reader.read(&mut record);
            if result.is_none() {
                break;
            }
            let _ = result.unwrap()?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if prev_name.as_ref() != Some(&qname) {
                prev_name.replace(qname);
                names_seen += 1;
            }
            if names_seen > n {
                self.read = Some(record.clone());
                break;
            }
            batch.push(record.clone());
        }

        Ok(batch)
    }

    /// Reads records up to and including the final record in a group of records with `name` as the qname.
    /// It's expected that the records are collated by qname, as will be the case for the bams written by
    /// this program to the workspace directory.
    pub fn take_until_qname(&mut self, name: &str) -> Result<Vec<bam::Record>> {
        let mut batch = Vec::new();
        let mut record = bam::Record::new();
        let mut found_match = false;

        // Include any cached read if it exists
        if let Some(cached_record) = self.read.take() {
            let qname = std::str::from_utf8(cached_record.qname())?;
            batch.push(cached_record.clone()); // Collect the cached read regardless of name
            if qname == name {
                found_match = true;
            }
        }

        // Process remaining records
        while let Some(result) = self.reader.read(&mut record) {
            let _ = result?;
            let qname = std::str::from_utf8(record.qname())?;

            // Check if the match condition has been met
            if qname == name {
                found_match = true;
                batch.push(record.clone());
            } else if found_match {
                // Stop collecting and cache the first mismatching read after a match
                self.read = Some(record.clone());
                break;
            } else {
                // Continue collecting reads until we find a match
                batch.push(record.clone());
            }
        }

        if !found_match {
            return Err(anyhow!("Target qname '{}' not found", name));
        }

        Ok(batch)
    }
}
