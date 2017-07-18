# bulk ATAC-seq pipeline

```
usage: bulk_ATAC-seq_pipeline.py [-h] [-p1 PAIRED1] [-p2 PAIRED2] -o OUTPUT
                                 [-n NAME] [-t THREADS] [-m MEMORY]
                                 [-q QUALITY] [-ref REFERENCE]
                                 [--picard_mark_dup PICARD_MARK_DUP]
                                 [--skip_trim] [--skip_align] [--skip_peaks]
                                 [--skip_bdg] [--skip_cleanup]

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
                        Maximum memory per thread for samtools sort [8]
  -q QUALITY, --quality QUALITY
                        Mapping quality cutoff for samtools [30]
  -ref REFERENCE, --reference REFERENCE
                        Path to reference genome
                        [/home/joshchiou/references/ucsc.hg19.fasta]
  --picard_mark_dup PICARD_MARK_DUP
                        Path to picard MarkDuplicates.jar
                        [/home/joshchiou/bin/MarkDuplicates.jar]

Skip processing steps:
  --skip_trim           Skip adapter trimming step
  --skip_align          Skip read alignment step
  --skip_peaks          Skip calling peaks step
  --skip_bdg            Skip making genome browser track
  --skip_cleanup        Skip cleanup operations
```
