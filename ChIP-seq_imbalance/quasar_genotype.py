#!/usr/bin/env python3
#===============================================================================
# quasar-genotype.py
#===============================================================================

"""Infer genotypes from ChIP-seq or ATAC-seq data using QuASAR"""




# Imports ----------------------------------------------------------------------

import argparse
import functools
import os
import os.path
import pickle
import subprocess
import sys
import socket
import tempfile

from multiprocessing import Pool

hostname = socket.gethostname()
if hostname == 'gatsby.ucsd.edu':
    sys.path.append('/home/data/kglab-python3-modules')
elif hostname == 'holden':
    sys.path.append('/lab/kglab-python3-modules')

import quasar
import seqalign
import wasp




# Function definitions ---------------------------------------------------------

def prepare_quasar_input(
    input_file_path: str,
    bam_dir: str,
    intermediate_dir: str,
    reference_genome_path: str,
    phred_quality_score: int,
    blacklist_path: str,
    snps_path: str,
    processes: int,
    memory: int,
    paired_end: bool,
    skip_preprocessing: bool,
    write_bam: bool,
    algorithm_switch_bp=70,
    algorithm=None
) -> str:
    """Format data into input files for QuASAR
    
    Parameters
    ----------
    input_file_path : str
        Path to an input file
    bam_dir : str
        Directory to write BAM files
    intermediate_dir : str
        Directory to write intermediate pileup / bed files
    reference_genome_path : str
        Path to reference genome
    phred_quality_score : int
        Minimum quality score for filtering alignment
    blacklist_path : str
        Path to ENCODE mappability blacklist
    snps_path : str
        Path to file containing SNPs to genotype
    processes : int
        Number of processes
    memory : int
        Memory limit
    paired_end : bool
        Indicator for paired_end reads
    skip_preprocessing : bool
        Indicator to skip preprocessing steps
    write_bam : bool
        Indicator to write a BAM file to disk
    algorithm_switch_bp : int
        Read length threshold for switching to `bwa mem`
    algorithm : str or None
        Force use of either `aln` or `mem` algorithm, if supplied
    
    Returns
    -------
    str
        Path to a QuASAR input file
    """
    input_file_path = (
        input_file_path.split(',')
        if
        ',' in input_file_path
        else
        input_file_path
    )
    bam_prefix, intermediate_prefix = (
        os.path.join(
            directory,
            (
                os.path.basename(
                    input_file_path
                    if
                    isinstance(input_file_path, str)
                    else
                    input_file_path[0]
                )
                .replace('.bam', '')
                .replace('.fastq', '')
                .replace('.gz', '')
            )
        )
        for
        directory
        in
        (bam_dir, intermediate_dir)
    )
    with open('{}.align.log'.format(bam_prefix), 'w') as log:
        sa = seqalign.SequenceAlignment(
            input_file=input_file_path,
            phred_quality_score=phred_quality_score,
            processes=processes,
            log=log,
            aligner=seqalign.BWA(
                reference_genome_path=reference_genome_path,
                trim_qual=15,
                algorithm_switch_bp=algorithm_switch_bp,
                algorithm=algorithm
            ),
            dedupper=wasp.RmDup(processes=processes, paired_end=paired_end)
        )
        if not skip_preprocessing:
            sa.apply_quality_filter()
            sa.remove_supplementary_alignments()
            sa.remove_blacklisted_reads(blacklist_path=blacklist_path)
            sa.samtools_sort(memory_limit=memory * processes)
            sa.samtools_index()
            sa.remove_mitochondrial_reads()
            sa.samtools_index()
            if write_bam:
                sa.write('{}.filt.bam'.format(bam_prefix))
            sa.remove_duplicates()
            sa.samtools_sort()
        compressed_pileup_bed_path = (
            '{}.pileup.bed.gz'.format(intermediate_prefix)
        )
        quasar.write_compressed_pileup_bed(
            sa.samtools_mpileup(
                positions=snps_path,
                reference_genome=reference_genome_path
            ),
            compressed_pileup_bed_path,
            snps_bed_path=snps_path
        )
    quasar.bed_to_quasar(compressed_pileup_bed_path)
    quasar_input_file_path = '{}.quasar.in.gz'.format(intermediate_prefix)
    return (
        quasar_input_file_path
        if
        os.path.isfile(quasar_input_file_path)
        else
        None
    )


def get_genotypes(
    input_list: list,
    bam_dir: str,
    intermediate_dir: str,
    reference_genome_path: str,
    phred_quality_score: int,
    blacklist_path: str,
    snps_path: str,
    processes: int,
    memory: int,
    paired_end: bool,
    skip_preprocessing: bool,
    write_bam: bool,
    algorithm_switch_bp=70,
    algorithm=None
):
    """Obtain genotypes from sequencing data using QuASAR
    
    Parameters
    ----------
    input_list : list
        List of input files
    bam_dir : str
        Directory to write BAM files
    intermediate_dir : str
        Directory to write intermediate pileup / bed files
    reference_genome_path : str
        Path to reference genome
    phred_quality_score : int
        Minimum quality score for filtering alignment
    blacklist_path : str
        Path to ENCODE mappability blacklist
    snps_path : str
        Path to file containing SNPs to genotype
    processes : int
        Number of processes
    memory : int
        Memory limit
    paired_end : bool
        Indicator for paired_end reads
    skip_preprocessing : bool
        Indicator to skip preprocessing steps
    write_bam : bool
        Indicator to write a BAM file to disk
    algorithm_switch_bp : int
        Read length threshold for switching to `bwa mem`
    algorithm : str or None
        Force use of either `aln` or `mem` algorithm, if supplied
    """
    
    n_input_files = len(input_list)
    with Pool(processes=min(processes, n_input_files)) as (
        pool
    ), tempfile.TemporaryDirectory(dir='/home/data/tmp') as (
        temp_dir_name
    ):
        return quasar.genotype(
            *filter(
                None,
                pool.map(
                    functools.partial(
                        prepare_quasar_input,
                        bam_dir=(
                            bam_dir
                            if
                            bam_dir
                            else
                            temp_dir_name
                        ),
                        intermediate_dir=(
                            intermediate_dir
                            if
                            intermediate_dir
                            else
                            temp_dir_name
                        ),
                        reference_genome_path=reference_genome_path,
                        phred_quality_score=phred_quality_score,
                        blacklist_path=blacklist_path,
                        snps_path=snps_path,
                        processes=max(1, int(processes / n_input_files)),
                        memory=memory,
                        paired_end=paired_end,
                        skip_preprocessing=skip_preprocessing,
                        write_bam=write_bam,
                        algorithm_switch_bp=algorithm_switch_bp,
                        algorithm=algorithm
                    ),
                    input_list
                )
            )
        )


def main(args):
    vcf = quasar.genotype_to_vcf(
        get_genotypes(
            args.input,
            args.bam_dir,
            args.inter_dir,
            args.reference,
            args.quality,
            args.blacklist,
            args.snps,
            args.processes,
            args.memory,
            args.paired_end,
            args.skip_preprocessing,
            args.write_bam,
            args.algorithm_switch_bp,
            args.algorithm
        ),
        sample_name=args.sample,
        snps_bed_path=args.snps,
        threshold=args.threshold,
        het_only=args.het_only
    )
    if args.vcf_chr and not args.quiet:
        vcf = tuple(vcf)
    if args.vcf_chr:
        quasar.write_split_vcf(vcf, args.vcf_chr)
    if not args.quiet:
        for line in vcf:
            print(line)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Infer genotypes from ChIP-seq or ATAC-seq data using QuASAR'
        )
    )
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument(
        'input',
        metavar='<path/to/sequencing_data.{fa/fq/bam}>',
        nargs='+',
        help='Paths to input FASTQ or BAM files'
    )
    io_group.add_argument(
        '--bam-dir',
        metavar='<path/to/bam_dir/>',
        help='directory in which to place BAM files'
    )
    io_group.add_argument(
        '--inter-dir',
        metavar='<path/to/inter_dir/>',
        help='prefix for intermediate files'
    )
    io_group.add_argument(
        '--vcf-chr',
        metavar='<output/vcf/prefix>',
        help='Prefix for output VCFs split by chromosome'
    )
    io_group.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress printing of VCF file to standard output'
    )
  
    align_group = parser.add_argument_group('Alignment arguments')
    align_group.add_argument(
        '--processes',
        metavar='<int>',
        type=int,
        default=4,
        help='Number of processes to use [4]'
    )
    align_group.add_argument(
        '--memory',
        metavar='<int>',
        type=int,
        default=8,
        help='Maximum memory per thread in GB [8]'
    )
    align_group.add_argument(
        '--quality',
        metavar='<int>',
        type=int,
        default=10,
        help='Mapping quality cutoff for samtools [10]'
    )
    align_group.add_argument(
        '--reference',
        metavar='<path/to/reference_genome.fa>',
        default='/home/joshchiou/references/ucsc.hg19.fasta',
        help=(
            'Path to reference genome prepared for BWA '
            '[/home/joshchiou/references/ucsc.hg19.fasta]'
        )
    )
    align_group.add_argument(
        '--blacklist',
        metavar='<path/to/blacklist.bed>',
        default='/home/data/encode/ENCODE.hg19.blacklist.bed',
        help=(
            'Path to ENCODE blacklist file '
            '[/home/data/encode/ENCODE.hg19.blacklist.bed]'
        )
    )
    align_group.add_argument(
        '--paired-end',
        action='store_true',
        help='Use for paired-end reads'
    )
    align_group.add_argument(
        '--write-bam',
        action='store_true',
        help='Write bam files to disk'
    )
    align_group.add_argument(
        '--algorithm-switch-bp',
        metavar='<int>',
        default=70,
        help='Read length threshold for switching to `bwa mem` [70]'
    )
    align_group.add_argument(
        '--algorithm',
        choices={'aln', 'mem'},
        default=None,
        help='Force use of either `bwa aln` or bwa mem`'
    )
    
    quasar_group = parser.add_argument_group('QuASAR arguments')
    quasar_group.add_argument(
        '--snps',
        metavar='<path/to/snps_file.bed>',
        default='/home/data/QuASAR/1KG_SNPs_filt.bed',
        help=(
            'BED file containing 1KGP SNPs '
            '[/home/data/QuASAR/1KG_SNPs_filt.bed]'
        )
    )
    quasar_group.add_argument(
        '--skip-preprocessing',
        action='store_true',
        help='skip preprocessing steps'
    )
    
    vcf_group = parser.add_argument_group('VCF arguments')
    vcf_group.add_argument(
        '--sample',
        metavar='<sample_id>',
        default='SAMPLE',
        help='Name for the sample [SAMPLE]'
    )
    vcf_group.add_argument(
        '--threshold',
        metavar='<float>',
        default=0.99,
        type=float,
        help='Probability threshold for genotype calls [0.99]'
    )
    vcf_group.add_argument(
        '--het-only',
        action='store_true',
        help='Output heterozygous variants only'
    )
    return parser.parse_args()




# Execute ----------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
