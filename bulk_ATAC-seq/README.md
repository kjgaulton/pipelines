# bulk ATAC-seq pipeline

NOTE: files necessary for the pipeline can be downloaded from the ./etc directory.

```
usage: bulk_ATAC-seq_pipeline.py [-h] [-p1 PAIRED1] [-p2 PAIRED2] -o OUTPUT
                                 [-n NAME] [-t THREADS] [-m MEMORY]
                                 [-q QUALITY] [-ref REFERENCE]
                                 [--picard_mark_dup PICARD_MARK_DUP]
                                 [--tss TSS] [--blacklist BLACKLIST]
                                 [--skip_trim] [--skip_align] [--skip_peaks]
                                 [--skip_bdg] [--skip_qc] [--skip_cleanup]

Pipeline for ATAC-seq to trim reads, align them to a reference genome, call
peaks, and generate a genome browser signal track.

optional arguments:
  -h, --help            show this help message and exit

I/O arguments:
  -p1 PAIRED1, --paired1 PAIRED1
                        Path to paired reads (1)
  -p2 PAIRED2, --paired2 PAIRED2
                        Path to paired reads (2)
  -o OUTPUT, --output OUTPUT
                        Output directory for processed files
  -n NAME, --name NAME  Output sample name to prepend

Alignment arguments:
  -t THREADS, --threads THREADS
                        Number of threads to use [4]
  -m MEMORY, --memory MEMORY
                        Maximum memory (in Gb) per thread for samtools sort
                        [8]
  -q QUALITY, --quality QUALITY
                        Mapping quality cutoff for samtools [30]
  -ref REFERENCE, --reference REFERENCE
                        Path to reference genome
                        [/home/joshchiou/references/ucsc.hg19.fasta]
  --picard_mark_dup PICARD_MARK_DUP
                        Path to picard MarkDuplicates.jar
                        [/home/joshchiou/bin/MarkDuplicates.jar]

QC arguments:
  --tss TSS             Path to TSS definitions for calculating ATAC signal
                        enrichment around TSS [/home/joshchiou/references/hg19
                        _gencode_tss_unique.bed]
  --blacklist BLACKLIST
                        Path to blacklist BED file to ignore ENCODE high
                        signal regions
                        [/home/joshchiou/references/ENCODE.hg19.blacklist.bed]

Skip processing steps:
  --skip_trim           Skip adapter trimming step [OFF]
  --skip_align          Skip read alignment step [OFF]
  --skip_peaks          Skip calling peaks step [OFF]
  --skip_bdg            Skip making genome browser track [OFF]
  --skip_qc             Skip ATAC qc step using ataqv [OFF]
  --skip_cleanup        Skip cleanup operations [OFF]

```
