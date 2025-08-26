use std::{collections::VecDeque, env, path::Path, process::Command};

use crate::{
    read_types::MappedReadPair,
    utils::{FifoGuard, ProcessGuard},
};
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
            result.unwrap()?;
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
            result?;
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

pub struct CollatedBamReader {
    reader: bam::Reader,
    current_batch: Option<MappedReadPair>,
    buffer: VecDeque<bam::Record>,
    _fifo_guard: FifoGuard,
    process_guard: ProcessGuard,
    exhausted: bool,
}

impl CollatedBamReader {
    pub fn new(bam_path: &Path, temp_dir: Option<&Path>, threads: Option<usize>) -> Result<Self> {
        // Get appropriate temp directory
        let temp_dir = match temp_dir {
            Some(dir) => dir.to_path_buf(),
            None => env::temp_dir(),
        };

        // Generate a unique FIFO name
        let pid = std::process::id();
        let timestamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();
        let fifo_name = format!("collated_bam_{}_{}.fifo", timestamp, pid);
        let fifo_path = temp_dir.join(fifo_name);

        // Create the FIFO with automatic cleanup
        let fifo_guard = FifoGuard::new(&fifo_path)?;

        // Set up threading
        let threads = match threads {
            Some(0) => 1,
            Some(t) => t,
            None => 1,
        };

        // Start the samtools process
        let process = Command::new("samtools")
            .arg("collate")
            .arg("-u") // uncompressed BAM
            .arg("-f") // Primary alignments only
            .arg("-@")
            .arg(threads.to_string())
            .arg("-o")
            .arg(&fifo_path)
            .arg(bam_path)
            .spawn()?;

        let process_guard = ProcessGuard::new(process);

        // Open the BAM file
        let reader = bam::Reader::from_path(&fifo_path)?;

        Ok(CollatedBamReader {
            reader,
            current_batch: None,
            buffer: VecDeque::new(),
            _fifo_guard: fifo_guard,
            process_guard,
            exhausted: false,
        })
    }

    // Helper method to fill the buffer with reads
    fn fill_buffer(&mut self) -> Result<bool> {
        let mut record = bam::Record::new();
        match self.reader.read(&mut record) {
            Some(Ok(())) => {
                self.buffer.push_back(record);
                Ok(true)
            }
            Some(Err(e)) => Err(anyhow!("Error reading BAM record: {}", e)),
            None => {
                // No more records - check if process completed successfully
                if !self.exhausted {
                    self.exhausted = true;
                    // Wait for process to complete and check status
                    self.process_guard.wait()?;
                }
                Ok(false)
            }
        }
    }

    // Process reads with the same name
    fn process_next_batch(&mut self) -> Result<Option<MappedReadPair>> {
        // If buffer is empty, try to fill it
        if self.buffer.is_empty() && !self.fill_buffer()? {
            return Ok(None); // No more records
        }

        // Get the first record to determine the read name
        let first_record = self.buffer.pop_front().unwrap();
        let read_name = String::from_utf8_lossy(first_record.qname()).to_string();

        // Create a new read pair
        let mut read_pair = MappedReadPair::new(&read_name);
        read_pair.insert(first_record)?;

        // Process all reads with the same name
        loop {
            // Fill buffer if needed
            if self.buffer.is_empty() && !self.fill_buffer()? {
                break; // No more records to process
            }

            // Check if next record has the same name
            if self.buffer.is_empty() {
                break;
            }

            let peek_name = String::from_utf8_lossy(self.buffer[0].qname()).to_string();
            if peek_name != read_name {
                break; // Different read name, end of batch
            }

            // Add record to the batch
            let record = self.buffer.pop_front().unwrap();
            read_pair.insert(record)?;
        }

        Ok(Some(read_pair))
    }

    // Check if the samtools process has completed
    pub fn process_finished(&mut self) -> Result<bool> {
        self.process_guard.try_wait()
    }
}

impl Iterator for CollatedBamReader {
    type Item = Result<MappedReadPair>;

    fn next(&mut self) -> Option<Self::Item> {
        // Check if we've already encountered an error with the subprocess
        if !self.exhausted && self.process_guard.try_wait().is_err() {
            return Some(Err(anyhow!("samtools process failed")));
        }

        match self.process_next_batch() {
            Ok(Some(pair)) => Some(Ok(pair)),
            Ok(None) => {
                // Double-check process status when we run out of records
                if !self.exhausted {
                    self.exhausted = true;
                    if let Err(e) = self.process_guard.wait() {
                        return Some(Err(e));
                    }
                }
                None
            }
            Err(e) => Some(Err(e)),
        }
    }
}
