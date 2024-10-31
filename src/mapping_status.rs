use anyhow::{anyhow, Result};
use std::collections::HashSet;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use crate::utils::block_filter;

#[derive(Debug, Clone)]
pub enum MappingStatus {
    Mapped, // All positions map to the reference (indels + mismatches permitted, also split/supplementary reads are merged)
    Unmapped, // Read is designated unmapped by the mapper
    PartiallyUnmapped(HashSet<usize>), // Some positions map, while others are clipped
    Suspicious, // All positions are unmapped, but the mapper didn't designate it as unmapped
    Unknown, // Reads have not been aligned yet, or other ambiguous status
}

fn get_clipped_positions(record: &bam::Record) -> Vec<usize> {
    if record.is_unmapped() {
        return (0..record.seq_len()).collect();
    }

    let mut positions = Vec::new();
    let mut pos = 0; // keeps track of position in the read

    let cigar = record.cigar();
    let cigar_ops: Box<dyn Iterator<Item = &Cigar>> = if record.is_reverse() {
        Box::new(cigar.iter().rev())
    } else {
        Box::new(cigar.iter())
    };

    for cigar_op in cigar_ops {
        match cigar_op {
            Cigar::SoftClip(len) | Cigar::HardClip(len) => {
                positions.extend(pos..pos + *len as usize);
                pos += *len as usize;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                pos += *len as usize;
            }
            _ => {}
        }
    }
    positions
}


pub fn get_mapping_status(records: &[bam::Record]) -> Result<MappingStatus> {
    if records.is_empty() {
        return Err(anyhow!("No records to process"));
    }

    let read_name = records[0].qname();
    let first_or_second = records[0].is_first_in_template();
    let len = records[0].seq_len();

    if records.iter().all(|record| record.is_unmapped()) {
        return Ok(MappingStatus::Unmapped);
    }

    let mut positions: HashSet<usize> =
        get_clipped_positions(&records[0]).iter().cloned().collect();

    for record in &records[1..] {
        if record.qname() != read_name {
            return Err(anyhow!("Read names do not match"));
        }

        if record.is_first_in_template() != first_or_second {
            return Err(anyhow!("Reads are not all the same position in the pair"));
        }

        if record.seq_len() != len {
            return Err(anyhow!("Reads are not all the same length"));
        }

        let p: HashSet<usize> = get_clipped_positions(record).iter().cloned().collect();

        if p.is_empty() {
            return Ok(MappingStatus::Mapped);
        }

        positions = positions.intersection(&p).cloned().collect();

        if positions.is_empty() {
            return Ok(MappingStatus::Mapped);
        }
    }
    if positions.is_empty() {
        return Ok(MappingStatus::Mapped);
    }

    if positions.len() == len {
        return Ok(MappingStatus::Suspicious);
    }

    // Mark read as partially unmapped if it passes block size and average block quality checks
    let positions = block_filter(positions, &records[0].qual(), 10, 20.0); // arbitrary values for min block size and quality

    if positions.is_empty() {
        return Ok(MappingStatus::Mapped);
    }

    Ok(MappingStatus::PartiallyUnmapped(positions))
}

