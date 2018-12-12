# Pipelines for analyzing snATAC-seq data

A typical experiment can be analyzed by using scripts/notebooks in the order:  
* demultiplex.py
* snATAC_pipeline.py
* snATAC_analysis.ipynb

demultplex.py - demultiplex barcodes and prepend them to the read name
snATAC_pipeline.py - process demultiplexed data and create count matrices using windows genome-wide
snATAC_analysis.ipynb - use count matrices to cluster cells from one experiment or across experiments

```
usage: snATAC_pipeline.py [-h] -r1 READ1 -r2 READ2 [-o OUTPUT] -n NAME
                          [-t THREADS] [-m MEMORY] [-q MAP_QUALITY]
                          [-ref REFERENCE] [--picard PICARD]
                          [--extension EXTENSION]
                          [--minimum-reads MINIMUM_READS]
                          [--window-size WINDOW_SIZE]
                          [--chrom-sizes CHROM_SIZES] [--blacklist BLACKLIST]
                          [--skip-trim] [--skip-align] [--skip-rmdup]
                          [--skip-matrix]

Align demulitplexed snATAC-seq reads to a reference genome.

optional arguments:
  -h, --help            show this help message and exit

I/O arguments:
  -r1 READ1, --read1 READ1
                        Paired-end reads file 1
  -r2 READ2, --read2 READ2
                        Paired-end reads file 2
  -o OUTPUT, --output OUTPUT
                        Output directory to store processed files
  -n NAME, --name NAME  Prefix for naming all output files

Alignment arguments:
  -t THREADS, --threads THREADS
                        Number of threads to use for alignment [8]
  -m MEMORY, --memory MEMORY
                        Maximum amount of memory (G) per thread for samtools
                        sort [4]
  -q MAP_QUALITY, --map_quality MAP_QUALITY
                        Mapping quality score filter for samtools [30]
  -ref REFERENCE, --reference REFERENCE
                        Path to the BWA-prepared reference genome

Remove duplicates arguments:
  --picard PICARD       Path to picard.jar

Matrix generation arguments:
  --extension EXTENSION
                        Read extension length
  --minimum-reads MINIMUM_READS
                        Minimum number of reads for barcode inclusion
  --window-size WINDOW_SIZE
                        Size (kb) to use for defining windows of accessibility
  --chrom-sizes CHROM_SIZES
                        Chromosome sizes file from UCSC
  --blacklist BLACKLIST
                        BED file of blacklisted regions

Skip steps:
  --skip-trim           Skip adapter trimming step
  --skip-align          Skip read alignment step
  --skip-rmdup          Skip duplicate removal step
  --skip-matrix         Skip matrix generation step
```
