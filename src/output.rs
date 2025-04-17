use crate::utils::{create_bgzf_fastq_writer, FastqWriter};
use anyhow::{anyhow, Result};
use bio::io::fastq;
use std::collections::HashMap;

pub struct OutputWriter {
    base_dir: std::path::PathBuf,
    files: HashMap<String, FastqWriter>,
    written_count: usize,
    fragment_count: usize,
}

impl OutputWriter {
    pub fn new(base_dir: std::path::PathBuf) -> Result<Self> {
        let mut this = Self {
            base_dir,
            files: HashMap::new(),
            written_count: 0,
            fragment_count: 0,
        };
        if this.base_dir.exists() {
            return Err(anyhow!(
                "Output directory already exists: {}",
                this.base_dir.display()
            ));
        }
        std::fs::create_dir_all(&this.base_dir)?;
        if !this.base_dir.is_dir() {
            return Err(anyhow!(
                "Problem creating output directory: {}",
                this.base_dir.display()
            ));
        }
        this.create_subdirectories()?;
        this.create_output_files()?;
        Ok(this)
    }

    fn create_subdirectories(&self) -> Result<()> {
        let subdirs = vec![
            "unmapped_pairs",
            "mapped_unmapped_pairs",
            "mapped_partial_pairs",
            "unmapped_partial_pairs",
            "partial_pairs",
            "unmapped_fragments",
        ];
        for subdir in subdirs {
            let path = self.base_dir.join(subdir);
            let result = std::fs::create_dir(&path);
            if let Err(e) = result {
                return Err(anyhow!(
                    "Problem creating subdirectory: {}: {}",
                    path.display(),
                    e
                ));
            }
        }
        Ok(())
    }

    fn create_output_files(&mut self) -> Result<()> {
        self.files.insert(
            "uuu1".to_string(),
            create_bgzf_fastq_writer(&self.base_dir.join("unmapped_pairs/u1.fq.gz"))?,
        );
        self.files.insert(
            "uuu2".to_string(),
            create_bgzf_fastq_writer(&self.base_dir.join("unmapped_pairs/u2.fq.gz"))?,
        );
        self.files.insert(
            "mum1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_unmapped_pairs/mapped_unmapped_m1.fq.gz"),
            )?,
        );
        self.files.insert(
            "muu2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_unmapped_pairs/mapped_unmapped_u2.fq.gz"),
            )?,
        );
        self.files.insert(
            "muu1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_unmapped_pairs/unmapped_mapped_u1.fq.gz"),
            )?,
        );
        self.files.insert(
            "mum2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_unmapped_pairs/unmapped_mapped_m2.fq.gz"),
            )?,
        );
        self.files.insert(
            "mpm1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_partial_pairs/mapped_partial_m1.fq.gz"),
            )?,
        );
        self.files.insert(
            "mpp2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_partial_pairs/mapped_partial_p2.fq.gz"),
            )?,
        );
        self.files.insert(
            "mpp1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_partial_pairs/partial_mapped_p1.fq.gz"),
            )?,
        );
        self.files.insert(
            "mpm2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("mapped_partial_pairs/partial_mapped_m2.fq.gz"),
            )?,
        );
        self.files.insert(
            "upu1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("unmapped_partial_pairs/unmapped_partial_u1.fq.gz"),
            )?,
        );
        self.files.insert(
            "upp2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("unmapped_partial_pairs/unmapped_partial_p2.fq.gz"),
            )?,
        );
        self.files.insert(
            "upp1".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("unmapped_partial_pairs/partial_unmapped_p1.fq.gz"),
            )?,
        );
        self.files.insert(
            "upu2".to_string(),
            create_bgzf_fastq_writer(
                &self
                    .base_dir
                    .join("unmapped_partial_pairs/partial_unmapped_u2.fq.gz"),
            )?,
        );
        self.files.insert(
            "ppp1".to_string(),
            create_bgzf_fastq_writer(&self.base_dir.join("partial_pairs/p1.fq.gz"))?,
        );
        self.files.insert(
            "ppp2".to_string(),
            create_bgzf_fastq_writer(&self.base_dir.join("partial_pairs/p2.fq.gz"))?,
        );
        self.files.insert(
            "frag".to_string(),
            create_bgzf_fastq_writer(&self.base_dir.join("unmapped_fragments/frag.fq.gz"))?,
        );
        Ok(())
    }

    pub fn write(&mut self, key: &str, record: &fastq::Record) -> Result<()> {
        let writer = self
            .files
            .get_mut(key)
            .ok_or_else(|| anyhow!("Invalid key: {}", key))?;
        writer.write_record(record)?;
        self.written_count += 1;
        Ok(())
    }

    pub fn write_batch(&mut self, key: &str, records: &[fastq::Record]) -> Result<()> {
        let writer = self
            .files
            .get_mut(key)
            .ok_or_else(|| anyhow!("Invalid key: {}", key))?;
        for record in records {
            writer.write_record(record)?;
            self.fragment_count += 1;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        for writer in self.files.values_mut() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn written(&self) -> usize {
        self.written_count
    }

    pub fn fragments(&self) -> usize {
        self.fragment_count
    }
}
