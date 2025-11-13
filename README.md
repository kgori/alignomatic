## Alignomatic
alignomatic is a program to extract reads from short read sequencing data that are unmapped, or partially unmapped, relative to a set of reference genomes.

The idea is to start with a pair of fastq files, for paired read1 and read2 reads, or better yet, a prealigned BAM file, and progressively align them to
a set of 'filter' reference genomes.

### Requirements

- If starting with fastq files, they need to be collated so that paired reads are at the same position in each file.
- The reference files need to be indexed with `bwa index`.

### Program options:

Options can be passed to the program through command line arguments, or via a JSON config file. If you use command line arguments, then these will be saved
into a JSON file for reuse / for your records / for fun.

alignomatic writes checkpoints of its progress, so if a run fails there's a chance you can pick it up from the point of failure and avoid recomputing everything
from the start.

```
Usage: alignomatic [OPTIONS]

Options:
      --config <CONFIG>
          Optional path to a config file (JSON format)
  -1, --fastq-first <FILE>

  -2, --fastq-second <FILE>

  -b, --bam-input <FILE>
          Optional pre-aligned BAM file input. Overrides fastq_first and fastq_second if provided.
  -i, --index <FILES>
          Reference files to map against. Must be Fasta format, and must have set of BWA index files.
  -o, --output-folder <FILE>
          Output folder for all fastq files. Will be created if it doesn't exist.
  -n, --batch-size <BATCH_SIZE>
          Batch size to process, in base pairs per thread. Default is 10000000.
  -t, --threads <THREADS>
          Number of threads to use. One thread will be reserved for the main program; any extra threads will be used for read mapping. Default is 1.
      --min-block-size <MIN_BLOCK_SIZE>
          Minimum size of a block of bases that will be considered unmapped. Default is 30.
      --min-block-quality <MIN_BLOCK_QUALITY>
          Minimum average base quality of a block of bases that will be considered unmapped. Default is 10.
  -h, --help
          Print help
```

### Output
alignomatic produces 16 fastq output files. These all have names like "reads_uu.1.fq.gz". The naming system is as follows:

`reads_XY.N.fa.gz`
  
  - files with the same `reads_XY` pattern contain paired reads.
  - `N` is the pairing information - 1 means the reads are first in pair, and 2 means they are second in pair.
  - `X` identifies the type of read that read 1 of the pair is.
  - `Y` identifies the type of read that read 2 of the pair is.
  - `X` and `Y` take one of the following codes:
     - `m` = mapped. This read maps fully to at least one reference genome
     - `u` = unmapped. This read is totally unmapped to any reference genome.
     - `f` = fragmentary mapping. This read maps partially to at least one reference genome. It contains a contiguous block of bases that don't map anywhere. The size and average base quality of this block are tunable through the program options.

`reads_uu.*.fq.gz` are totally foreign DNA, as far as the filter references go. 
The other files' reads are potentially useful for anchoring insertions of foreign DNA into a reference genome.

#### Name
alignomatic is a stupid name, but I watched too much Rocko's Modern Life growing up. All machines are "-o-matic" to me.
