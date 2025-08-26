use crate::cli::ProgramOptions;
use crate::mapping_status::get_mapping_status;
use crate::mapping_status::MappingStatus::*;
use crate::read_pair_io::ReadPairIterator;
use crate::read_types::MappedReadPair;
use anyhow::{anyhow, Result};
use bio::io::fastq;
use log::{debug, warn};
use rust_htslib::bam;
use std::collections::HashSet;
use std::fs::File;
use std::os::unix::io::AsRawFd;
use std::path::{Component, Path, PathBuf};
use std::process::Command;

/// Intercepts any stderr output from the wrapped function
pub fn silence_stderr<T, F>(f: F) -> Result<T>
where
    F: FnOnce() -> T,
{
    let dev_null = File::open("/dev/null")?;
    let dev_null_fd = dev_null.as_raw_fd();
    let stderr_fd = unsafe { libc::dup(libc::STDERR_FILENO) };

    if stderr_fd < 0 {
        return Err(anyhow!("Failed to duplicate stderr file descriptor"));
    }

    if unsafe { libc::dup2(dev_null_fd, libc::STDERR_FILENO) } < 0 {
        return Err(anyhow!("Failed to redirect stderr to /dev/null"));
    }

    let result = f();

    if unsafe { libc::dup2(stderr_fd, libc::STDERR_FILENO) } < 0 {
        return Err(anyhow!("Failed to restore stderr"));
    }

    Ok(result)
}

#[derive(Debug)]
struct Block {
    /// Represents a contiguous range of positions. The range is inclusive on the start and exclusive on the end.
    start: usize,
    end: usize,
}

impl Block {
    fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    fn len(&self) -> usize {
        self.end - self.start
    }
}

/// Converts a hashset of usize positions into sorted vector of contiguous blocks
fn to_blocks(set: HashSet<usize>) -> Vec<Block> {
    if set.is_empty() {
        return Vec::new();
    }

    let mut sorted = set.into_iter().collect::<Vec<usize>>();
    sorted.sort();

    let mut blocks = Vec::<Block>::new();
    let mut start: usize = sorted[0];
    let mut curr: usize = start;

    for next in sorted[1..].iter() {
        if *next - curr == 1 {
            curr = *next;
        } else {
            blocks.push(Block::new(start, curr + 1));
            start = *next;
            curr = start;
        }
    }
    blocks.push(Block::new(start, curr + 1));

    blocks
}

/// View the unmapped positions as contiguous blocks, and accept a block as valid unmapped positions
/// if it passes length and average quality requirements
pub fn block_filter(
    positions: HashSet<usize>,
    qual: &[u8],
    min_block_size: usize,
    min_avg_qual: f32,
) -> HashSet<usize> {
    let blocks = to_blocks(positions);
    let mut filtered_positions: HashSet<usize> = HashSet::new();
    for block in blocks {
        if block.len() < min_block_size {
            continue;
        }

        let block_qual = (block.start..block.end)
            .map(|pos| qual[pos] as f32)
            .sum::<f32>()
            / block.len() as f32;
        if block_qual >= min_avg_qual {
            filtered_positions.extend(block.start..block.end);
        }
    }
    filtered_positions
}

pub trait Slice {
    fn slice(&self, range: std::ops::Range<usize>) -> Result<Self>
    where
        Self: Sized;
}

impl Slice for fastq::Record {
    fn slice(&self, range: std::ops::Range<usize>) -> Result<Self> {
        if range.end > self.seq().len() {
            return Err(anyhow!("Range exceeds sequence length"));
        }
        // Annotate the record description with the slice info
        let desc = match self.desc() {
            Some(d) => format!("{}_{}-{}", d, range.start, range.end),
            None => format!("{}-{}", range.start, range.end),
        };
        Ok(fastq::Record::with_attrs(
            self.id(),
            Some(&desc),
            &self.seq()[range.clone()],
            &self.qual()[range],
        ))
    }
}

pub fn fastq_to_unmapped_fragments(
    record: &fastq::Record,
    positions: &HashSet<usize>,
) -> Result<Vec<fastq::Record>> {
    let blocks = to_blocks(positions.clone());
    let mut fragments = Vec::new();
    for block in blocks {
        fragments.push(record.slice(block.start..block.end)?);
    }
    Ok(fragments)
}

pub fn check_file_exists(file: &Path) -> Result<()> {
    let existence = file.try_exists();
    match existence {
        Ok(true) => Ok(()),
        _ => Err(anyhow!("File not found: {}", file.to_string_lossy())),
    }
}

pub fn check_directory_exists(dir: &Path) -> Result<()> {
    let existence = dir.try_exists();
    match (existence, dir.is_dir()) {
        (Ok(true), true) => Ok(()),
        (Ok(true), false) => Err(anyhow!(
            "Path is not a directory: {}",
            dir.to_string_lossy()
        )),
        _ => Err(anyhow!("Directory not found: {}", dir.to_string_lossy())),
    }
}

pub fn bam_to_fastq(record: &bam::Record) -> Result<fastq::Record> {
    if record.is_secondary() || record.is_supplementary() {
        return Err(anyhow!("Secondary or supplementary alignment"));
    }
    let seq = record.seq().as_bytes();
    let qual = record.qual().iter().map(|q| q + 33).collect::<Vec<u8>>();
    let id = String::from_utf8(record.qname().to_vec())?;

    if record.is_reverse() {
        let seq = seq
            .iter()
            .map(|base| bio::alphabets::dna::complement(*base))
            .rev()
            .collect::<Vec<u8>>();
        let qual = qual.into_iter().rev().collect::<Vec<u8>>();
        return Ok(fastq::Record::with_attrs(&id, None, &seq, &qual));
    }
    Ok(fastq::Record::with_attrs(&id, None, &seq, &qual))
}

pub fn normalize_path<P: AsRef<Path>>(path: P) -> Result<PathBuf> {
    let mut stack = Vec::new();
    let path = std::env::current_dir()?.join(path);

    for component in path.components() {
        match component {
            Component::ParentDir => {
                stack.pop();
            }
            Component::CurDir => {}
            _ => stack.push(component),
        }
    }

    let mut normalized = PathBuf::new();
    for component in stack {
        normalized.push(component);
    }

    Ok(normalized)
}

pub fn estimate_results_batch_size(opts: &ProgramOptions) -> Result<usize> {
    debug!(target: "UTILS", "Estimating batch size");
    let fastq1 = opts.fastq_first.clone().unwrap();
    let fastq2 = opts.fastq_second.clone().unwrap();
    let threads = opts.threads;
    let batch_size = opts.batch_size;
    let bases = threads * batch_size;
    let mut read_pair_iter = ReadPairIterator::new(fastq1, fastq2)?;
    debug!(target: "UTILS", "Attempting to load {} bases", bases);
    let read_pairs = read_pair_iter.take_bases(bases);
    let n_read_pairs = read_pairs.len();
    debug!(target: "UTILS", "Loaded {} read pairs", n_read_pairs);
    if n_read_pairs == 0 {
        Err(anyhow!("No read pairs found in the input files"))
    } else {
        Ok(n_read_pairs)
    }
}

pub fn read_pair_is_unmapped(mapped_pair: &MappedReadPair, opts: &ProgramOptions) -> Result<bool> {
    let read1_status = get_mapping_status(&mapped_pair.read1, opts)?;
    let read2_status = get_mapping_status(&mapped_pair.read2, opts)?;
    match (read1_status, read2_status) {
        (Mapped, Mapped) => Ok(false),
        (Suspicious, _) | (_, Suspicious) => {
            warn!("Read pair {} has a suspicious alignment", mapped_pair.id());
            Ok(false)
        }
        _ => Ok(true),
    }
}

// RAII struct to ensure cleanup of FIFO on drop
pub struct FifoGuard {
    path: PathBuf,
}

impl FifoGuard {
    pub fn new(path: &Path) -> Result<Self> {
        if path.exists() {
            std::fs::remove_file(path)?;
        }

        let status = Command::new("mkfifo")
            .arg(path)
            .status()
            .map_err(|e| anyhow!("Failed to create FIFO: {}", e))?;

        if !status.success() {
            return Err(anyhow!("mkfifo failed with status: {}", status));
        }

        Ok(FifoGuard {
            path: path.to_path_buf(),
        })
    }
}

impl Drop for FifoGuard {
    fn drop(&mut self) {
        if self.path.exists() {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

pub struct ProcessGuard {
    process: std::process::Child,
    finished: bool,
}

impl ProcessGuard {
    pub fn new(process: std::process::Child) -> Self {
        ProcessGuard {
            process,
            finished: false,
        }
    }

    pub fn try_wait(&mut self) -> Result<bool> {
        if self.finished {
            return Ok(true);
        }

        match self.process.try_wait() {
            Ok(Some(status)) => {
                self.finished = true;
                if status.success() {
                    Ok(true)
                } else {
                    Err(anyhow!("Process exited with error: {}", status))
                }
            }
            Ok(None) => Ok(false),
            Err(e) => Err(anyhow!("Failed to wait for process: {}", e)),
        }
    }

    pub fn wait(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }

        let status = self.process.wait()?;

        self.finished = true;

        if !status.success() {
            Err(anyhow!("Process exited with error: {}", status))
        } else {
            Ok(())
        }
    }
}

impl Drop for ProcessGuard {
    fn drop(&mut self) {
        if !self.finished {
            // Try to wait for the process to finish gracefully
            if self.process.wait().is_err() {
                // If waiting fails, kill the process
                let _ = self.process.kill();
            }
        }
    }
}

pub fn extract_primary_read(records: &[bam::Record]) -> Result<fastq::Record> {
    if records.is_empty() {
        return Err(anyhow!("No records found"));
    }
    let primary_read = records
        .iter()
        .find(|record| !record.is_secondary() && !record.is_supplementary())
        .map(bam_to_fastq);

    match primary_read {
        Some(Ok(read)) => Ok(read),
        Some(Err(e)) => Err(anyhow!("Error converting BAM record to FASTQ: {}", e)),
        None => {
            for record in records {
                eprintln!("Record: {:?}", record);
            }
            Err(anyhow!("No primary alignment found"))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::read_pair_io::ReadPairIterator;
    use rust_htslib::bam::Read;

    #[test]
    fn test_bam_to_fastq() {
        let mut bam_reader = bam::Reader::from_path("data/test.bam").unwrap();
        let mut fastq_reader = ReadPairIterator::new(
            PathBuf::from("data/test.1.fq.gz"),
            PathBuf::from("data/test.2.fq.gz"),
        )
        .unwrap();

        let mut converted_1 = Vec::new();
        let mut converted_2 = Vec::new();
        let reference = fastq_reader.take_pairs(10000);

        loop {
            let bam_record = bam_reader.records().next().unwrap().unwrap();
            if bam_record.is_secondary() || bam_record.is_supplementary() {
                continue;
            }
            let first = bam_record.is_first_in_template();
            let last = bam_record.is_last_in_template();
            let converted = bam_to_fastq(&bam_record).unwrap();
            if first && converted_1.len() < reference.len() {
                converted_1.push(converted.clone());
            }
            if last && converted_2.len() < reference.len() {
                converted_2.push(converted);
            }

            if converted_1.len() == reference.len() && converted_2.len() == reference.len() {
                break;
            }
        }

        for (i, (r1, r2)) in converted_1.iter().zip(converted_2.iter()).enumerate() {
            assert_eq!(r1.id(), reference[i].id);
            assert_eq!(r2.id(), reference[i].id);
            assert_eq!(r1.seq(), reference[i].read1.seq());
            assert_eq!(r2.seq(), reference[i].read2.seq());
            assert_eq!(r1.qual(), reference[i].read1.qual());
            assert_eq!(r2.qual(), reference[i].read2.qual());
        }
    }
}
