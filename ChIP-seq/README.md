# ChIP-seq pipeline
```
usage: ChIP-seq_pipeline.py [-h] -t TREATMENT -c CONTROL -o OUTPUT [-n NAME]
                            [-p PROCESSES] [-m MEMORY] [-q QUALITY] -ref
                            REFERENCE -markdup MARKDUP [--broad]
                            [--skip_align] [--skip_peaks] [--skip_track]

Pipeline for ChIP to align reads to a reference genome, and then call peaks.

optional arguments:
  -h, --help            show this help message and exit

I/O arguments:
  -t TREATMENT, --treatment TREATMENT
                        Path to treatment .fastq.gz
  -c CONTROL, --control CONTROL
                        Path to control .fastq.gz
  -o OUTPUT, --output OUTPUT
                        Output directory for processed files
  -n NAME, --name NAME  Output sample name to prepend

Alignment and rmdup arguments:
  -p PROCESSES, --processes PROCESSES
                        Number of processes to use
  -m MEMORY, --memory MEMORY
                        Maximum memory per thread
  -q QUALITY, --quality QUALITY
                        Mapping quality cutoff for samtools
  -ref REFERENCE, --reference REFERENCE
                        Path to reference genome prepared for BWA
  -markdup MARKDUP, --markdup MARKDUP
                        Path to MarkDuplicates.jar

MACS2 parameters:
  --broad               Broad peak option for MACS2 callpeak

Skip processing:
  --skip_align          Skip read alignment step
  --skip_peaks          Skip calling peaks step
  --skip_track          Skip making signal track for genome browser
```
