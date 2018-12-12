#!/usr/bin/env python3
#===============================================================================
# wasp_map.py
#===============================================================================

"""Implementation of the re-mapping procedure to eliminate reference bias 
detailed in WASP.

See: https://github.com/bmvdgeijn/WASP/tree/master/mapping
"""



# Imports ----------------------------------------------------------------------

import argparse
import functools
import math
import os
import os.path
import socket
import sys

from multiprocessing import Pool

hostname = socket.gethostname()
if hostname == 'gatsby.ucsd.edu':
    sys.path.append('/home/data/kglab-python3-modules')
elif hostname == 'holden':
    sys.path.append('/lab/kglab-python3-modules')

import seqalign
import wasp




# Function definitions ---------------------------------------------------------

def setup_directory(output_dir, snp_dir):
    """Create directories for the WASP mapping pipeline as needed
    
    Parameters
    ----------
    output_dir : str
        path to output directory for the pipeline
    snp_dir : str
        path to directory where SNP files will be stored
    """
    
    for directory in (output_dir, snp_dir) + tuple(
        os.path.join(output_dir, subdir)
        for
        subdir
        in
        (
            'map1',
            'find_intersecting_snps',
            'map2',
            'filter_remapped_reads',
            'pileup'
        )
    ):
        if not os.path.isdir(directory):
            os.mkdir(directory)


def map_reads(file_path, trim_qual, processes=1, algorithm=None):
    """Map some reads
    
    Parameters
    ----------
    file_path : str, list, tuple
        Path to reads on disk that will be mapped. A string for single-end
        reads, a list or tuple containing two strings for paired-end reads.
    trim_qual : int
        MAPQ score for BWA read trimming
    processes : int
        Maximum number of processes
    
    Returns
    -------
    seqalign.SequenceAlignment
        Object representing the aligned data
    """
    
    if len(file_path) == 2:
        if file_path[1] == 'paired':
            file_path = file_path[0]
    with seqalign.SequenceAlignment(
        file_path,
        processes=processes,
        aligner=seqalign.BWA(
            trim_qual=trim_qual,
            algorithm=algorithm,
            algorithm_switch_bp=100
        )
    ) as sa:
        sa.samtools_sort()
        sa.samtools_index
        return sa


def preprocess_sample(
    sample_name,
    file_path,
    output_dir,
    processes=1,
    algorithm=None
):
    """Apply preprocessing steps to an input sample
    
    Writes a preprocessed BAM file to the map1 subdir of the output directory.
    If the input file path matches the output file path, it is assumed that
    preprocessing has already been done and this step is skipped.
    
    Parameters
    ----------
    sample_name : str
        Sample name
    file_path : list
        List containin one or two paths to input files (for single- or
        paired-end reads)
    output_dir : str
        Output directory for the pipeline
    processes : int
        Maximum number of processes to use for this sample
    """
    
    map1_path = os.path.join(
        output_dir,
        'map1',
        '{}.sort.bam'.format(sample_name)
    )
    file_path = file_path[0] if len(file_path) == 1 else file_path
    
    if file_path != map1_path:   
        alignment = map_reads(
            file_path,
            trim_qual=15,
            processes=processes,
            algorithm=algorithm
        )
        alignment.write(map1_path)


def preprocessing_step(sample_list, output_dir, processes=1, algorithm=None):
    """Apply preprocessing steps to all input samples
    
    Preprocesses multiple samples in parallel
    
    Parameters
    ----------
    sample_list : list
        List of sample information from the argument parser
    output_dir : str 
        Output directory for the pipeline
    processes : int
        Maximum number of processes to use
    """
    
    n_samples = len(sample_list)
    with Pool(processes=min(processes, n_samples)) as pool:
        pool.starmap(
            functools.partial(
                preprocess_sample,
                output_dir=output_dir,
                processes=max(1, math.floor(processes / n_samples)),
                algorithm=algorithm
            ),
            (
                (sample_name, input_file_path)
                for
                sample_name, *input_file_path
                in
                sample_list
            )
        )


def _find_intersecting_snps(
    input_file_path,
    bam_filename,
    snp_dir,
    output_dir=None
):
    """Run WASP's `find_intersecting_snps.py` script
    
    This function wraps wasp.find_intersecting_snps and handles endedness by
    interpreting the type of the provided input_file_path
    
    Parameters
    ----------
    input_file_path : str, list, tuple
        Path to an input file (or paths to two files if paired-end) as provided
        to the argument parser. The length of this parameter is checked to
        determine endedness.
    bam_filename : str, bytes
        Path to an input BAM file (str), or an input BAM file in memory (bytes)
    snp_dir : str
        Path to directory containing SNP files
    output_dir : str
        Path to directory where output files will be written
    """
    
    wasp.find_intersecting_snps(
        bam_filename,
        snp_dir,
        is_paired_end=len(input_file_path) == 2,
        is_sorted=True,
        output_dir=output_dir
    )


def find_intersecting_snps(
    sample_list,
    output_dir,
    snp_dir,
    processes=1,
    memory_limit=5
):
    """Run WASP's `find_intersecting_snps.py` script
    
    This function runs _find_intersecting_snps on multiple files in parallel
    
    Parameters
    ----------
    sample_list : list
        List of sample information from the argument parser
    output_dir : str 
        Output directory for the pipeline
    snp_dir : str
        Path to directory containing SNP files
    processes : int
        Maximum number of processes to use
    memory_limit : int
        Approximate memory limit in gigabytes [5]
    """
    
    with Pool(
        processes=min(
            processes,
            math.floor(memory_limit / 4),
            len(sample_list)
        )
    ) as pool:
        pool.starmap(
            functools.partial(
                _find_intersecting_snps,
                snp_dir=snp_dir,
                output_dir=os.path.join(
                    output_dir,
                    'find_intersecting_snps'
                )
            ),
            (
                (
                    input_file_path,
                    os.path.join(
                        output_dir,
                        'map1',
                        '{}.sort.bam'.format(sample_name)
                    )
                )
                for
                sample_name, *input_file_path
                in
                sample_list
            )
        )


def remap_sample(
    sample_name,
    file_path,
    output_dir,
    processes=1,
    algorithm=None
):
    """Remap an input sample
    
    Writes a remapped BAM file to the map2 subdirectory of the output directory
    
    Parameters
    ----------
    sample_name : str
        Sample name
    file_path : list
        List containin one or two paths to input files (for single- or
        paired-end reads)
    output_dir : str
        Output directory for the pipeline
    processes : int
        Maximum number of processes to use for this sample
    """
    
    alignment = map_reads(
        file_path,
        trim_qual=0,
        processes=processes,
        algorithm=algorithm
    )
    alignment.write(
        os.path.join(output_dir, 'map2', '{}.sort.bam'.format(sample_name))
    )


def remapping_step(sample_list, output_dir, processes=1, algorithm=None):
    """Remap all input samples
    
    Remaps multiple samples in parallel
    
    Parameters
    ----------
    sample_list : list
        List of sample information from the argument parser
    output_dir : str 
        Output directory for the pipeline
    processes : int
        Maximum number of processes to use
    """
    
    n_samples = len(sample_list)
    with Pool(processes=min(processes, n_samples)) as pool:
        pool.starmap(
            functools.partial(
                remap_sample,
                output_dir=output_dir,
                processes=max(1, math.floor(processes / n_samples)),
                algorithm=algorithm
            ),    
            (
                (
                    sample_name,
                    (
                        tuple(
                            os.path.join(
                                output_dir,
                                'find_intersecting_snps',
                                '{}.sort.remap.fq{}.gz'.format(
                                    sample_name,
                                    end
                                )
                            )
                            for
                            end
                            in
                            (1, 2)
                        )
                        if
                        len(input_file_path) == 2
                        else
                        os.path.join(
                            output_dir,
                            'find_intersecting_snps',
                            '{}.sort.remap.fq.gz'.format(sample_name)
                        )
                    )
                )
                for
                sample_name, *input_file_path
                in
                sample_list
            )
        )


def filter_remapped_reads(sample_list, output_dir, processes=1, memory_limit=5):
    """Run WASP's `filter_remapped_reads.py` script
    
    This function runs wasp.filter_remapped_reads on multiple files in parallel
    
    Parameters
    ----------
    sample_list : list
        List of sample information from the argument parser
    output_dir : str 
        Output directory for the pipeline
    processes : int
        Maximum number of processes to use
    memory_limit : int
        Approximate memory limit in gigabytes [5]
    """
    
    with Pool(
        processes=min(processes, memory_limit, len(sample_list))    
    ) as pool:
        pool.starmap(
            wasp.filter_remapped_reads,
            (
                (
                    os.path.join(
                        output_dir,
                        'find_intersecting_snps',
                        '{}.sort.to.remap.bam'.format(sample_name)
                    ),
                    os.path.join(
                        output_dir,
                        'map2',
                        '{}.sort.bam'.format(sample_name)
                    ),
                    os.path.join(
                        output_dir,
                        'filter_remapped_reads',
                        '{}.keep.bam'.format(sample_name)
                    )
                )
                for
                sample_name, *input_file_path
                in
                sample_list
            )
        )


def merge_and_rmdup(*sequence_alignments, paired_end=False, processes=1):
    """Merge kept and remapped reads and apply WASP's dedupper
    
    Parameters
    ----------
    *sequence_alignments
        One or more seqalign.SequenceAlignment objects
    paired_end : bool
        If True, `rmdup_pe.py` will be used, otherwise `rmdup.py` will be used
    processes : int
        Maximum number of processes for samtools sort
    """
    
    sa = seqalign.merge(
        *sequence_alignments,
        processes=processes,
        dedupper=wasp.RmDup(paired_end=paired_end, processes=processes)
    )
    sa.samtools_sort()
    sa.samtools_index()
    sa.remove_supplementary_alignments()
    sa.remove_duplicates()
    sa.samtools_sort()
    sa.samtools_index()
    return sa


def merge_rmdup_pileup(
    sample_name,
    input_file_path,
    output_dir,
    snp_dir,
    reference_genome_path,
    processes=1
):
    """Merge kept and remapped reads, remove duplicates, and generate a pileup
    
    The pileup will be written to the pileup subdir of the output directory
    
    Parameters
    ----------
    sample_name : str
        Sample name
    input_file_path : str, list, tuple
        Path to an input file (or paths to two files if paired-end) as provided
        to the argument parser. The length of this parameter is checked to
        determine endedness.
    output_dir : str 
        Output directory for the pipeline
    snp_dir : str
        Path to directory containing SNP files
    reference_genome_path : str
        Path to a reference genome on disk
    """
    
    alignment = merge_and_rmdup(
        os.path.join(
            output_dir,
            'filter_remapped_reads',
            '{}.keep.bam'.format(sample_name)
        ),
        os.path.join(
            output_dir,
            'find_intersecting_snps',
            '{}.sort.keep.bam'.format(sample_name)
        ),
        paired_end=len(input_file_path) == 2,
        processes=processes
    )
    with open(
        os.path.join(
            output_dir,
            'pileup',
            '{}.pileup'.format(sample_name)
        ),
        'w'
    ) as f:
        f.write(
            alignment.samtools_mpileup(
                positions=os.path.join(snp_dir, 'snps.positions.txt'),
                reference_genome=reference_genome_path
            )
            .decode()
        )


def merge_rmdup_pileup_steps(
    sample_list,
    output_dir,
    snp_dir,
    reference_genome_path,
    processes=1
):
    """Merge reads, remove duplicates, and generate a pileup for all samples
    
    Pileups will be generated for multiple samples in parallel
    
    Parameters
    ----------
    sample_list : list
        List of sample information from the argument parser
    output_dir : str 
        Output directory for the pipeline
    snp_dir : str
        Path to directory containing SNP files
    reference_genome_path : str
        Path to a reference genome on disk
    """
    
    n_samples = len(sample_list)
    with Pool(processes=min(processes, n_samples)) as pool:
        pool.starmap(
            functools.partial(
                merge_rmdup_pileup,
                output_dir=output_dir,
                snp_dir=snp_dir,
                reference_genome_path=reference_genome_path,
                processes=max(1, math.floor(processes / n_samples))
            ),
            (
                (sample_name, input_file_path)
                for
                sample_name, *input_file_path
                in
                sample_list
            )
        )


def remove_intermediate_files(samples, output_dir):
    """Remove intermediate files generated during the procedure
    
    Parameters
    ----------
    samples : list
        The list of samples obtained from command line arguments
    output_dir : str
        Directory for output files
    """
    
    for file_path in (
        file_path
        for
        sample_name
        in
        (sample[0] for sample in samples)
        for
        file_path
        in
        (
            os.path.join(
                output_dir,
                'map2',
                '{}.sort.bam'.format(sample_name)
            ),
            os.path.join(
                output_dir,
                'find_intersecting_snps',
                '{}.sort.remap.fq.gz'.format(sample_name)
            ),
            os.path.join(
                output_dir,
                'find_intersecting_snps',
                '{}.sort.to.remap.bam'.format(sample_name)
            ),
            os.path.join(
                output_dir,
                'find_intersecting_snps',
                '{}.sort.keep.bam'.format(sample_name)
            ),
            os.path.join(
                output_dir,
                'filter_remapped_reads',
                '{}.keep.bam'.format(sample_name)
            )
        ) + tuple(
            os.path.join(
                output_dir,
                'find_intersecting_snps',
                '{}.sort.remap.fq{}.gz'.format(
                    sample_name,
                    end
                )
            )
            for
            end
            in
            (1, 2)
        )
        if
        os.path.isfile(file_path)
    ):
        os.remove(file_path)


def remove_empty_directories(output_dir):
    """Remove empty directories left over after the pipeline is completed
    
    Parameters
    ----------
    output_dir : str
        Directory for output files
    """
    
    for directory in (
        os.path.join(output_dir, subdir)
        for
        subdir
        in
        ('find_intersecting_snps', 'map2', 'filter_remapped_reads')
        if
        os.path.isdir(os.path.join(output_dir, subdir))
        if
        not os.listdir(os.path.join(output_dir, subdir))
    ):
        os.rmdir(directory)

def main(args):
    """The main pipeline"""
    
    # Step 0: set up the output directory
    setup_directory(args.output_dir, args.snp_dir)
    
    # Step 1: create SNP files
    if args.vcf_format and args.vcf_sample:
        wasp.write_positions_snps(
            args.vcf_format,
            (
                os.path.join(args.snp_dir, 'snps')
                if
                args.snp_dir
                else
                os.path.join(args.output_dir, 'snp_dir', 'snps')
            ),
            het_only=True,
            r2=args.r2,
            samples=args.vcf_sample,
            keep_filtered_vcfs=True
        )
    
    # Stop here if no sample was given
    if args.sample:
    
        # Step 2: preprocess sequencing data
        preprocessing_step(
            args.sample,
            args.output_dir,
            processes=args.processes,
            algorithm=args.algorithm
        )
    
        # Step 3: find intersecting snps
        find_intersecting_snps(
            args.sample,
            args.output_dir,
            args.snp_dir,
            processes=args.processes,
            memory_limit=args.memory_limit
        )
    
        # Step 4: remap
        remapping_step(
            args.sample,
            args.output_dir,
            processes=args.processes,
            algorithm=args.algorithm
        )
    
        # Step 5: filter remapped reads
        filter_remapped_reads(
            args.sample,
            args.output_dir,
            processes=args.processes,
            memory_limit=args.memory_limit
        )
    
        # Step 6-7+: Merge remapped_reads, remove duplicates, generate pileups
        merge_rmdup_pileup_steps(
            args.sample,
            args.output_dir,
            args.snp_dir,
            args.reference_genome,
            processes=args.processes
        )
    
        # Clean up intermediate files
        if not args.save_intermediate:
            remove_intermediate_files(args.sample, args.output_dir)
        
    # Clean up empty directories
    if not args.save_intermediate:
        remove_empty_directories(args.output_dir)


def main_skip_to_6(args):
    # Stop here if no sample was given
    if args.sample:
    
        # Step 6-7+: Merge remapped_reads, remove duplicates, generate pileups
        merge_rmdup_pileup_steps(
            args.sample,
            args.output_dir,
            args.snp_dir,
            args.reference_genome,
            processes=args.processes
        )
    
        # Clean up intermediate files
        if not args.save_intermediate:
            remove_intermediate_files(args.sample, args.output_dir)
        
    # Clean up empty directories
    if not args.save_intermediate:
        remove_empty_directories(args.output_dir)


# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Implementation of the re-mapping procedure to eliminate reference '
            'bias detailed in WASP.'
        )
    )
    parser.add_argument(
        'output_dir',
        metavar='<path/to/output_dir/>',
        nargs='?',
        default='.',
        help=(
            'Directory for output files. If this is not provided, the current '
            'directory is used.'
        )
    )
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument(
        '--sample',
        metavar=('<sample_name>', '<path/to/sequencing_data.{fa/fq/bam}>'),
        action='append',
        nargs='+',
        help=(
            'Sample name followed by file path (or paths if paired-end reads). '
            'If the input file is a BAM file with paired-end reads, add a '
            'third argument "paired".'
        )
    )
    io_group.add_argument(
        '--save-intermediate',
        action='store_true',
        help='Do not delete intermediate files'
    )
    config_group = parser.add_argument_group('config arguments')
    config_group.add_argument(
        '--vcf-format',
        metavar='<vcf/format/string.chr{}.vcf.gz>',
        help=(
            'Format for input VCF files. "{}" will be replaced by a '
            'chromosome number or character.'
        )
    )
    config_group.add_argument(
        '--vcf-sample',
        metavar='<sample_id>',
        help='Sample to use from VCF files'
    )
    config_group.add_argument(
        '--snp-dir',
        metavar='<path/to/snp_dir/>',
        help='Path to SNP directory'
    )
    config_group.add_argument(
        '--reference-genome',
        metavar='<path/to/reference_genome.fa>',
        default='/home/joshchiou/references/ucsc.hg19.fasta',
        help=(
            'Path to reference genome '
            '[/home/joshchiou/references/ucsc.hg19.fasta]'
        )
    )
    config_group.add_argument(
        '--r2',
        metavar='<float>',
        default=0.9,
        type=float,
        help=(
            'R2 imputation quality threshold [0.9]. Set this to 0 for '
            'non-imputed VCF data.'
        )
    )
    config_group.add_argument(
        '--algorithm',
        default=None,
        choices={'aln', 'mem'},
        help='Force BWA algorithm'
    )
    resource_group = parser.add_argument_group('resource arguments')
    resource_group.add_argument(
        '--processes',
        metavar='<int>',
        default=1,
        type=int,
        help='Maximum number of processes allowed [1]'
    )
    resource_group.add_argument(
        '--memory-limit',
        metavar='<int>',
        default=5,
        type=int,
        help='Approximate memory limit in gigabytes [5]'
    )
    troubleshooting_group = parser.add_argument_group(
        'troubleshooting arguments'
    )
    troubleshooting_group.add_argument(
        '--skip-to-6',
        action='store_true',
        help='Skip to step 6, merging the reads after remapping'
    )
    args = parser.parse_args()
    if not args.snp_dir:
        args.snp_dir = os.path.join(args.output_dir, 'snp_dir')
    return args




# Execute ----------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    if args.memory_limit < 5:
        raise Exception('Please provide at least 5 GB of memory')
    if not args.skip_to_6:
        main(args)
    else:
        main_skip_to_6(args)
