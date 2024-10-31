use anyhow::{anyhow, Result};
use bio::io::fastq;
use std::fs::File;
use std::collections::HashSet;
use std::os::unix::io::AsRawFd;

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

/// Conversion from rust_htslib::bam::Record to bio::io::fastq::Record
pub fn bam_to_fastq(
    record: &rust_htslib::bam::Record,
) -> Result<fastq::Record, Box<dyn std::error::Error>> {
    let name = std::str::from_utf8(record.qname())?;
    let desc = None;
    let seq = record.seq().as_bytes();
    let qual = record.qual().iter().map(|q| q + 33).collect::<Vec<u8>>();
    Ok(fastq::Record::with_attrs(name, desc, &seq, &qual))
}

/// Calculates the average base quality of a record
pub fn avg_quality(record: &rust_htslib::bam::Record) -> f32 {
    record.qual().iter().map(|x| *x as f32).sum::<f32>() / record.qual().len() as f32
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

    let mut sorted = set.into_iter()
        .collect::<Vec<usize>>();
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
pub fn block_filter(positions: HashSet<usize>, qual: &[u8], min_block_size: usize, min_avg_qual: f32) -> HashSet<usize> {
    let blocks = to_blocks(positions);
    let mut filtered_positions: HashSet<usize> = HashSet::new();
    for block in blocks {
        if block.len() < min_block_size {
            continue;
        }

        let block_qual = (block.start..block.end).map(|pos| qual[pos] as f32).sum::<f32>() / block.len() as f32;
        if block_qual >= min_avg_qual {
            filtered_positions.extend(block.start..block.end);
        }
    }
    filtered_positions
}

pub trait Slice {
    fn slice(&self, range: std::ops::Range<usize>) -> Result<Self> where Self: Sized;
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
            &self.qual()[range]))
    }
}

pub fn fastq_to_unmapped_fragments(record: &fastq::Record, positions: &HashSet<usize>) -> Result<Vec<fastq::Record>> {
    let blocks = to_blocks(positions.clone());
    let mut fragments = Vec::new();
    for block in blocks {
        fragments.push(record.slice(block.start..block.end)?);
    }
    Ok(fragments)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bam_to_fastq() {
        let mut bam_rec = rust_htslib::bam::Record::new();
        bam_rec.set(b"test", None, b"ACGT", &[0, 0, 0, 41]);

        let fastq_rec = fastq::Record::with_attrs("test", None, b"ACGT", b"!!!J");

        assert_eq!(bam_to_fastq(&bam_rec).unwrap(), fastq_rec,);
    }

    #[test]
    fn test_avg_quality() {
        let mut bam_rec = rust_htslib::bam::Record::new();
        bam_rec.set(b"test", None, b"ACGT", &[10, 12, 15, 45]);
        assert_eq!(avg_quality(&bam_rec), 20.5);
    }
}
