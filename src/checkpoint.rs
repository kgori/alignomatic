use anyhow::Result;
use log::debug;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

#[derive(Serialize, Deserialize, Debug)]
pub struct Checkpoint {
    reference: String,
    fastq_input_1: String,
    fastq_input_2: String,
    bam_output: String,
    fastq_output_1: String,
    fastq_output_2: String,
    batch_size: usize,
    min_block_size: usize,
    min_block_quality: f32,
}

fn blake3_hash(path: &PathBuf) -> Result<String> {
    debug!(target: "Checkpoint", "Hashing {:?}", path);
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut hasher = blake3::Hasher::new();
    let mut buffer = [0u8; 8192];

    loop {
        let n = reader.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }

    debug!(target: "Checkpoint", "Hashed {:?}", path);
    Ok(hasher.finalize().to_hex().to_string())
}

impl Checkpoint {
    pub fn create(
        reference: &PathBuf,
        fastq_input_1: &PathBuf,
        fastq_input_2: &PathBuf,
        bam_output: &PathBuf,
        fastq_output_1: &PathBuf,
        fastq_output_2: &PathBuf,
        batch_size: usize,
        min_block_size: usize,
        min_block_quality: f32,
    ) -> Result<Checkpoint> {
        let paths = vec![
            reference,
            fastq_input_1,
            fastq_input_2,
            bam_output,
            fastq_output_1,
            fastq_output_2,
        ];

        let threads = paths
            .into_iter()
            .map(|path| {
                let path = path.clone();
                std::thread::spawn(move || blake3_hash(&path))
            })
            .collect::<Vec<_>>();

        let hashes = threads
            .into_iter()
            .map(|thread| thread.join().unwrap())
            .collect::<Result<Vec<_>>>()?;

        Ok(Checkpoint {
            reference: hashes[0].clone(),
            fastq_input_1: hashes[1].clone(),
            fastq_input_2: hashes[2].clone(),
            bam_output: hashes[3].clone(),
            fastq_output_1: hashes[4].clone(),
            fastq_output_2: hashes[5].clone(),
            batch_size,
            min_block_size,
            min_block_quality,
        })
    }

    pub fn matches(&self, other: &Checkpoint) -> bool {
        self.reference == other.reference
            && self.fastq_input_1 == other.fastq_input_1
            && self.fastq_input_2 == other.fastq_input_2
            && self.bam_output == other.bam_output
            && self.fastq_output_1 == other.fastq_output_1
            && self.fastq_output_2 == other.fastq_output_2
            && self.batch_size == other.batch_size
            && self.min_block_size == other.min_block_size
            && self.min_block_quality == other.min_block_quality
    }
}

pub fn write_checkpoint(checkpoint: &Checkpoint, path: &PathBuf) -> Result<()> {
    debug!(target: "Checkpoint", "Writing checkpoint to {:?}", path);
    let file = File::create(path)?;
    let writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(writer, checkpoint)?;
    Ok(())
}

pub fn read_checkpoint(path: &PathBuf) -> Result<Checkpoint> {
    debug!(target: "Checkpoint", "Reading checkpoint from {:?}", path);
    let file = File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let checkpoint: Checkpoint = serde_json::from_reader(reader)?;
    Ok(checkpoint)
}
