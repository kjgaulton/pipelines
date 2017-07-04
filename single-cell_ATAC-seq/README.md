# sc-ATAC
single cell ATAC-seq pipeline

## Introduction
**sc-ATAC** is an analysis pipeline for single cell ATAC-seq data based on the paper by [Cusanovich et al., 2015](http://science.sciencemag.org/content/348/6237/910).  
sc-ATAC can either be run as a pipeline or as independent steps.

```
sc-ATAC (single cell ATAC-seq pipeline by the Gaulton lab at UC San Diego)
Last Update:	06/22/2017
Maintained by:	Josh Chiou <joshchiou@ucsd.edu>
Gaulton Lab:	https://gaultonlab.org

Usage: sc-ATAC.pipeline.py [-h] [-i1 INDEX1] [-i2 INDEX2] [-p1 PAIRED1] [-p2 PAIRED2] [-o OUTPUT]
				[-n NAME] [-a ADAPTERS] [-t THREADS] [-mem MEMORY] [-ref REFERENCE]
				[-mark_dup PICARD_MARK_DUP] [--ambiguous_bases AMBIGUOUS_BASES]
				[--mismatches MISMATCHES] [--closest_distance CLOSEST_DISTANCE]
				[--min_reads MIN_READS] [--map_quality MAP_QUALITY]
				

Arguments:
  -h			--help			show this help message and exit.
  -i1			INDEX1			Path to index1 fastq file.
  -i2			INDEX2			Path to index2 fastq file.
  -p1			PAIRED1			Path to paired-end reads 1 fastq file.
  -p2			PAIRED2			Path to paired-end reads 2 fastq file.
  -o			OUTPUT			Output directory to store processed files.
  -n			NAME			Prefix for naming all output files.
  -a			ADAPTERS		Path to directory that contains Illumina adapters [../illumina_adapters/].
  -t			THREADS			Maximum amount of threads to be used [8].
  -mem			MEMORY			Maximum amount of memory to be used per thread [4].
  -ref			REFERENCE		Reference genome for aligning reads.
  -mark_dup		PICARD_MARK_DUP		Path to picard's MarkDuplicates.jar tool [../picard/MarkDuplicates.jar].
  --ambiguous_bases	AMBIGUOUS_BASES		Maximum number of ambiguous bases N per barcode [12].
  --mismatches		MISMATCHES		Maximum number of mismatches [3].
  --closest_distance	CLOSEST_DISTANCE	Minimum edit distance that secondary match needs to be away from the primary match [1].
  --min_reads		MIN_READS		Minimum number of reads for output barcodes [3].
  --map_quality		MAP_QUALITY		Mapping quality score filter for samtools [10].
```

## Installation
### Instructions
```
git clone https://github.com/joshchiou/sc-ATAC.git
cd sc-ATAC/
chmod u+x bin/*
```
### Dependencies
|Software	|Version	|Link									|
|:---		|:---		|:---									|
|Python 3	|3.4+		|https://www.python.org/downloads/release/python-361/			|
|Levenshtein	|0.12.0		|https://pypi.python.org/pypi/python-Levenshtein/			|
|trim_galore!	|0.4.4+		|https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/	|
|fastqc		|0.11.5		|https://www.bioinformatics.babraham.ac.uk/projects/fastqc/		|
|bwa		|0.7.12+	|https://github.com/lh3/bwa/releases/					|
|samtools	|1.12+		|http://www.htslib.org/download/					|

## Pipeline
The pipeline *sc-ATAC.pipeline.py* consists of 5 steps, each of which can be run individually.
1. *sc-ATAC.extract-bc.py*- 	extract barcodes from index files and adds the barcodes into the read names
2. *sc-ATAC.correct-bc.py* - 	correct barcodes for mismatches and filters reads based on barcode mismatches and edit distance
3. *sc-ATAC.align_reads.py* - 	trim reads using trim_galore and then aligns reads to the reference genome using bwa-mem
4. *sc-ATAC.split_cells.py* - 	split the aligned reads based on unique barcodes, then remove duplicates using Picard 
5. *sc-ATAC.make_matrix.py* - 	make a binary counts matrix in 200 bp windows

## Output
De-duplicated reads for each cell that passes the read filter.  
Summary for each step in the pipeline process.
