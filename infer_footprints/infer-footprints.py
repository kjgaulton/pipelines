#!/usr/bin/env python3

#------------------------------------------------------------------------------#
#                             infer-footprints.py                              #
#------------------------------------------------------------------------------#

# Script to streamline the peak-calling and footprinting pipeline




#--------------------------------- Imports ------------------------------------#

import argparse
import gzip
import os.path
import sys
import tempfile

sys.path.append('/home/data/kglab-python3-modules')
sys.path.append('/lab/kglab-python3-modules')

import chippeaks
import footprints
import seqalign
import hg19




#--------------------------- Function definitions -----------------------------#

def main(args):
    with seqalign.SequenceAlignment(
        input_file=(args.input, args.in2) if args.in2 else args.input,
        phred_quality_score=args.qual,
        processes=args.max_processes,
        aligner=seqalign.BwaAligner(
            reference_genome_path=args.reference_genome,
            trim=args.trim_reads
        )
    ) as sa:
        sa.remove_supplementary_alignments()
        sa.samtools_sort(memory_limit=args.memory_limit)
        sa.samtools_index()
        bam = sa.bam
        bam_index = sa.index
    
    with chippeaks.ChipPeaks(input_bam=bam) as cp:
        cp.blacklist(args.blacklist)
        peaks = cp.peaks
    
    with footprints.Footprints(
        input_bam=bam,
        input_bam_index=bam_index,
        input_peaks=peaks,
        reference_genome_path=args.reference_genome,
        motifs_db_path=args.motifs_db_path,
        centipede_path=args.centipede_path,
        processes=args.max_processes
    ) as fp:
        fp.cleans_up_footprints = False
        fp.dump_json(args.output)
        footprints_dict = fp.footprints




# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Pipeline for peak calling'
        )
    )
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument(
        '-i',
        '--input',
        required=True,
        help='Input file--reads or BAM'
    )
    required_group.add_argument(
        '-o',
        '--output',
        required=True,
        help='Output file--where to put the DNase footprints.'
    )
    required_group.add_argument(
        '--memory-limit',
        required=True,
        type=int,
        help='Approximate memory limit for samtools sort, in gigabytes'
    )
    required_group.add_argument(
        '--blacklist',
        required=True,
        help='Encode blacklist file'
    )
    required_group.add_argument(
        '--motifs-db-path',
        required=True,
        help='Path to motifs database'
    )
    optional_group = parser.add_argument_group('optional arguments')
    optional_group.add_argument(
        '-i2',
        '--in2',
        help='For paired-end reads'
    )
    optional_group.add_argument(
        '-q',
        '--qual',
        default='30',
        help='cutoff Phred score for filtering out FASTQ reads'
    )
    optional_group.add_argument(
        '--max-processes',
        type=int,
        default=1,
        help='Maximum number of processes allowed'
    )
    optional_group.add_argument(
        '--trim-reads',
        type=int,
        help='Trimming reads???? idek'
    )
    optional_group.add_argument(
        '--motifs-list',
        help='List of motifs'
    )
    optional_group.add_argument(
        '--centipede-path',
        help='Path to R script for centipede'
    )
    optional_group.add_argument(
        '--reference-genome',
        default=hg19.path,
        help='Path to reference genome'
    )

    args = parser.parse_args()

    return args




#--------------------------------- Execute ------------------------------------#

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
