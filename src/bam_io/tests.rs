use super::CollatedBamReader;
use anyhow::Result;
use std::path::PathBuf;

#[test]
fn test_collated_bam_reader_merge_logic() -> Result<()> {
    // Use the existing test BAM file
    let input_bam = PathBuf::from("data/test.sorted.bam");
    assert!(
        input_bam.exists(),
        "Test BAM file not found at data/test.sorted.bam"
    );

    // Use a small batch size to force chunking and merging, but not so small that we hit open file limits.
    // If the file has e.g. 20k records, batch_size=5 would create 4000 files, likely exceeding ulimit.
    // batch_size=2000 should be safe (creates ~10-50 files).
    let batch_size = 2000;

    let reader = CollatedBamReader::new(&input_bam, None, None, Some(batch_size))?;

    // Collect all records to verify order and integrity
    let mut collated_records = Vec::new();

    for result in reader {
        let mapped_pair = result?;

        // Check read pair structure
        assert!(
            !mapped_pair.read1.is_empty() || !mapped_pair.read2.is_empty(),
            "Read pair {} is empty",
            mapped_pair.id()
        );

        // Collect qnames to verify sorting
        if !mapped_pair.read1.is_empty() {
            collated_records.push(mapped_pair.read1[0].clone());
        } else if !mapped_pair.read2.is_empty() {
            collated_records.push(mapped_pair.read2[0].clone());
        }
    }

    // Verify we got some records
    assert!(!collated_records.is_empty(), "No records read from BAM");

    // Check sorting order (QNAME)
    // Note: CollatedBamReader returns MappedReadPairs which groups by name.
    // We just need to check that the names we see are non-decreasing.
    for i in 0..collated_records.len() - 1 {
        let qname_a = collated_records[i].qname();
        let qname_b = collated_records[i + 1].qname();
        assert!(
            qname_a <= qname_b,
            "Records not sorted: {:?} > {:?}",
            String::from_utf8_lossy(qname_a),
            String::from_utf8_lossy(qname_b)
        );
    }

    Ok(())
}
