# lung snATAC-seq pipeline

## Step 1: Run lung_snATAC_pipeline.py

```
usage: lung_snATAC_pipeline.py [-h] -r1 READ1 -r2 READ2 [-o OUTPUT] -n NAME
                               [-t THREADS] [-m MEMORY] [-q MAP_QUALITY]
                               [-ref REFERENCE] --picard PICARD
                               [--shift SHIFT] [--extsize EXTSIZE]
                               [--minimum-reads MINIMUM_READS]
                               [--window-size WINDOW_SIZE]
                               [--chrom-sizes CHROM_SIZES]
                               [--blacklist-file BLACKLIST_FILE]
                               [--promoter-file PROMOTER_FILE] [--skip-trim]
                               [--skip-align] [--skip-rmdup] [--skip-qc]
                               [--skip-matrix]

Process combinatorial barcoding snATAC-seq data.

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
                        Path to the BWA indexed reference genome

Remove duplicates arguments:
  --picard PICARD       Path to picard.jar

Matrix generation arguments:
  --shift SHIFT         Read shift size
  --extsize EXTSIZE     Read extension size
  --minimum-reads MINIMUM_READS
                        Minimum number of reads for barcode inclusion
  --window-size WINDOW_SIZE
                        Size (kb) to use for defining windows of accessibility
  --chrom-sizes CHROM_SIZES
                        Chromosome sizes file from UCSC
  --blacklist-file BLACKLIST_FILE
                        BED file of blacklisted regions
  --promoter-file PROMOTER_FILE
                        BED file of autosomal promoter regions

Skip steps:
  --skip-trim           Skip adapter trimming step
  --skip-align          Skip read alignment step
  --skip-rmdup          Skip duplicate removal step
  --skip-qc             Skip QC metrics calculation step
  --skip-matrix         Skip matrix generation step
```

## Step 2
Use the output files from Step 1 to cluster snATAC-seq data with the lung_snATAC.ipynb notebook.
