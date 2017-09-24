`python3 infer-footprints.py --help`

    usage: infer-footprints.py [-h] -i INPUT -o OUTPUT --memory-limit MEMORY_LIMIT
                               --blacklist BLACKLIST --motifs-db-path
                               MOTIFS_DB_PATH [-i2 IN2] [-q QUAL]
                               [--max-processes MAX_PROCESSES]
                               [--trim-reads TRIM_READS]
                               [--motifs-list MOTIFS_LIST]
                               [--centipede-path CENTIPEDE_PATH]
                               [--reference-genome REFERENCE_GENOME]

    Pipeline for peak calling

    optional arguments:
      -h, --help            show this help message and exit

    required arguments:
      -i INPUT, --input INPUT
                            Input file--reads or BAM
      -o OUTPUT, --output OUTPUT
                            Output file--where to put the DNase footprints.
      --memory-limit MEMORY_LIMIT
                            Approximate memory limit for samtools sort, in
                            gigabytes
      --blacklist BLACKLIST
                            Encode blacklist file
      --motifs-db-path MOTIFS_DB_PATH
                            Path to motifs database

    optional arguments:
      -i2 IN2, --in2 IN2    For paired-end reads
      -q QUAL, --qual QUAL  cutoff Phred score for filtering out FASTQ reads
      --max-processes MAX_PROCESSES
                            Maximum number of processes allowed
      --trim-reads TRIM_READS
                            Trimming reads???? idek
      --motifs-list MOTIFS_LIST
                            List of motifs
      --centipede-path CENTIPEDE_PATH
                            Path to R script for centipede
      --reference-genome REFERENCE_GENOME
                            Path to reference genome
