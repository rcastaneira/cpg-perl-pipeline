# cpg-perl-pipeline

A Perl pipeline for detecting CpG islands in genomic sequence and summarizing their methylation, built to reproduce and cross-check part of my comparative epigenomics work on the human chromosome 2 fusion locus.

## Current state

`cpg_islands.pl` detects CpG islands using the Takai and Jones (2002) criteria:

- length of at least 500 bp
- GC content of at least 55%
- observed/expected CpG ratio of at least 0.65

It slides a window along each sequence, merges overlapping windows that pass, re-checks each merged region as a whole, and writes BED-style output.

```
perl cpg_islands.pl --step 10 test.fa > islands.bed
```

Options: `--window` (default 500), `--step` (default 1), `--gc` (default 0.55), `--oe` (default 0.65).

`test.fa` contains an 800 bp CpG-rich block (positions 3000 to 3800) inside AT-rich background, plus a background-only sequence. Expected output is one island at 2770 to 3930 and nothing on the second sequence.

## Roadmap

1. Rolling-count window scan, so each step is constant time instead of recounting the whole window
2. FASTA reader as a reusable module (`lib/`) with `Test::More` tests
3. Run on the CHM13 chr2 fusion window (chr2:113,000,000 to 115,000,000) and compare with the island calls from my capstone analysis
4. Per-island methylation summary from modkit bedMethyl output
5. Documentation and example results

## Requirements

Perl 5 (core modules only so far: `strict`, `warnings`, `Getopt::Long`).
