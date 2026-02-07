use std::{
    collections::VecDeque,
    path::{Path, PathBuf},
    sync::mpsc::{sync_channel, Receiver, SyncSender},
    thread::{self, JoinHandle},
};

use crate::read_types::MappedReadPair;
use anyhow::{anyhow, Result};
use itertools::Itertools;
use log::info;
use rust_htslib::bam::{self, Read};
use tempfile::TempDir;

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
    receiver: Receiver<Result<bam::Record>>,
    current_batch: Option<MappedReadPair>,
    buffer: VecDeque<bam::Record>,
    join_handle: Option<JoinHandle<Result<()>>>,
    exhausted: bool,
}

impl CollatedBamReader {
    pub fn new(
        bam_path: &Path,
        _temp_dir: Option<&Path>,
        _threads: Option<usize>,
        batch_size: Option<usize>,
    ) -> Result<Self> {
        let bam_path_buf = bam_path.to_path_buf();
        let batch_size = batch_size.unwrap_or(1_000_000);

        // Create a bounded channel (sync_channel) to provide backpressure.
        // Size 5000 records seems reasonable?
        // 1M batch size in bytes vs records.
        // Let's use a channel bound of 2048 records. This prevents the sorter from reading
        // the entire file into memory if the consumer is slow.
        let (sender, receiver) = sync_channel(2048);

        // Spawn the collation thread
        let join_handle = thread::spawn(move || collate_bam_task(bam_path_buf, sender, batch_size));

        Ok(CollatedBamReader {
            receiver,
            current_batch: None,
            buffer: VecDeque::new(),
            join_handle: Some(join_handle),
            exhausted: false,
        })
    }

    // Helper method to fill the buffer with reads
    fn fill_buffer(&mut self) -> Result<bool> {
        match self.receiver.recv() {
            Ok(Ok(record)) => {
                self.buffer.push_back(record);
                Ok(true)
            }
            Ok(Err(e)) => Err(anyhow!("Error from collation thread: {}", e)),
            Err(_) => {
                // Channel closed (sender dropped), meaning task finished.
                if !self.exhausted {
                    self.exhausted = true;
                    // Check thread status for panic or error returned from the function
                    if let Some(handle) = self.join_handle.take() {
                        match handle.join() {
                            Ok(Ok(())) => {}
                            Ok(Err(e)) => return Err(anyhow!("Collation thread failed: {}", e)),
                            Err(_) => return Err(anyhow!("Collation thread panicked")),
                        }
                    }
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
}

impl Iterator for CollatedBamReader {
    type Item = Result<MappedReadPair>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.process_next_batch() {
            Ok(Some(pair)) => Some(Ok(pair)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

fn collate_bam_task(
    input_bam: PathBuf,
    sender: SyncSender<Result<bam::Record>>,
    batch_size: usize,
) -> Result<()> {
    info!(target: "IO", "Starting native BAM collation");

    // Create a temp dir for the sort chunks
    let temp_dir = TempDir::new()?;
    let mut chunk_files = Vec::new();

    // 1. Read, Sort, Chunk
    {
        let mut reader = bam::Reader::from_path(&input_bam)?;
        let header = bam::Header::from_template(reader.header());
        let mut records = Vec::with_capacity(batch_size);

        for record in reader.records() {
            let record = record?;
            records.push(record);

            if records.len() >= batch_size {
                let chunk_path = temp_dir
                    .path()
                    .join(format!("chunk_{}.bam", chunk_files.len()));
                sort_and_write_chunk(&mut records, &chunk_path, &header)?;
                chunk_files.push(chunk_path);
            }
        }

        // Flush remaining
        if !records.is_empty() {
            let chunk_path = temp_dir
                .path()
                .join(format!("chunk_{}.bam", chunk_files.len()));
            sort_and_write_chunk(&mut records, &chunk_path, &header)?;
            chunk_files.push(chunk_path);
        }
    }

    info!(target: "IO", "BAM collation: Merging {} chunks", chunk_files.len());

    if chunk_files.is_empty() {
        // Empty input, nothing to write
        return Ok(());
    }

    // 2. Merge

    // Create readers for all chunks
    // Note: bam::Reader is not Send, but we are in one thread.
    let mut chunk_readers = Vec::new();
    for path in &chunk_files {
        chunk_readers.push(bam::Reader::from_path(path)?);
    }

    // Merge logic
    // We use kmerge_by. We need to handle Result<Record>.
    // Comparison: if a.qname < b.qname.
    // If we encounter error during iteration, we panic or try to pass it through?
    // kmerge_by terminates on None.

    let iterators = chunk_readers.iter_mut().map(|r| r.records());

    let merged = iterators.kmerge_by(|res_a, res_b| {
        match (res_a, res_b) {
            (Ok(a), Ok(b)) => a.qname() < b.qname(),
            (Err(_), _) => true, // bubble error up? priority doesn't matter much if erroring
            (_, Err(_)) => false,
        }
    });

    for record_res in merged {
        let record = record_res?;
        if sender.send(Ok(record)).is_err() {
            // Receiver dropped, stop processing
            break;
        }
    }

    info!(target: "IO", "Native BAM collation complete");
    Ok(())
}

fn sort_and_write_chunk(
    records: &mut Vec<bam::Record>,
    path: &Path,
    header: &bam::Header,
) -> Result<()> {
    // Sort by QNAME
    records.sort_by(|a, b| a.qname().cmp(b.qname()));

    // Write
    // Use uncompressed BAM (Level 0) for speed in temp files?
    // rust-htslib Format::Bam defaults to compressed.
    // There isn't an easy way to set compression level in Writer::from_path without Config.
    // Let's just use standard BAM.
    let mut writer = bam::Writer::from_path(path, header, bam::Format::Bam)?;
    for record in records.drain(..) {
        writer.write(&record)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
