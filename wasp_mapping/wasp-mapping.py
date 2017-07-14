#!/usr/bin/env python3

#------------------------------------------------------------------------------#
#                             wasp-mapping.py                                  #
#------------------------------------------------------------------------------#

# Implementation of the re-mapping procedure to eliminate reference bias
# detailed in WASP.




#--------------------------------- Imports ------------------------------------#

# Import some necessary modules. We'll use subprocess to call things like bwa,
# anaconda-associated python, etc., Pool from multiprocessing to run
# certain things in parallel, and argparse to parse arguments.

import subprocess
from multiprocessing import Pool
import argparse




#--------------------------- Low-level functions ------------------------------#

# Run bwa aln with options for multithreading and read trimming. 
def bwa_aln(args, read_trimming, fastq, sai):
    with open(sai, 'w') as f:
        subprocess.call(
            [
                'bwa',
                'aln', 
                '-t', str(args.max_processes),
                '-q', str(read_trimming),
                args.reference_genome,
                fastq
            ],
            stdout=f
        )




# Run bwa samse
def bwa_samse(args, fastq, sai, sam):
    with open(sam, 'w') as f:
        subprocess.call(
            [
                'bwa',
                'samse',
                args.reference_genome,
                sai,
                fastq
            ],
            stdout=f
        )




# Run bwa sampe. The -P setting increases speed but requires 4-5 GB of memory
def bwa_sampe(args, fastq1, sai1, fastq2, sai2, sam):
    with open(sam, 'w') as f:
        subprocess.call(
            [
                'bwa',
                'sampe',
                '-P',
                args.reference_genome,
                sai1,
                sai2,
                fastq1,
                fastq2
            ],
            stdout=f
        )




# run STAR
def star(args, fastq_list, prefix):
    subprocess.call(
        [
            'STAR',
            '--runThreadN', str(args.max_processes),
            '--genomeDir', args.STAR_genomeDir,
            '--readFilesIn', fastq_list[0],
        ] + [
            fastq_list[-1]
        ] * (
            args.paired_end
        ) + [
            '--readFilesCommand', 'zcat',
            
            '--twopassMode', 'Basic',
            
            '--sjdbGTFfile', args.STAR_sjdbGTFfile,
            
            '--outFilterType', 'BySJout',
            '--outFilterMultimapNmax', '20',
            '--alignSJoverhangMin', '8',
            '--alignSJDBoverhangMin', '1',
            '--outFilterMismatchNmax', '999',
            '--alignIntronMin', '20',
            '--alignIntronMax', '1000000',
            '--alignMatesGapMax', '1000000',
            
            '--outSAMmapqUnique', '50',
            '--outMultimapperOrder', 'Random',
            '--outSAMmultNmax', '1',
            
            '--outSAMattributes', 'NH', 'HI', 'AS', 'nM', 'MD',
            
            '--outSAMunmapped', 'Within',
            '--outFilterMismatchNoverLmax', '0.04',
            '--sjdbScore', '1',
            '--genomeLoad', 'NoSharedMemory',
            '--outSAMheaderHD', '@HD VN:1.4 SO:unsorted',
            
            '--outFileNamePrefix', prefix
        ]
    )





# Run samtools view with BAM output and options for base quality filtering and
# multithreading
def samtools_view(args, sam, bam):
    with open(bam, 'w') as f:
        subprocess.call(
            [
                'samtools',
                'view',
                '-Sbq', '30',
                '-@', str(args.max_processes),
                sam
            ],
            stdout=f
        )




# Remove supplementary alignments from a BAM file
def remove_supplementary(args, bam, nosupbam):
    with open(nosupbam, 'w') as f:
        subprocess.call(
            [
                'samtools',
                'view',
                '-b',
                '-F', '0x800',
                '-@', str(args.max_processes),
                bam
            ],
            stdout=f
        )




# Run samtools sort with multithreading options
def samtools_sort(args, bam, sortedbam):
    subprocess.call(
        [
            'samtools',
            'sort',
            '-m', ''.join((str(int(1024 / 47 * args.memory_limit)), 'M')),
            '-@', str(args.max_processes),
            '-o', sortedbam,
            bam
        ]
    )




# Run samtools index
def samtools_index(bam):
    subprocess.call(['samtools', 'index', bam])




#--------------------------- Mid-level functions ------------------------------#

# Initial bwa aln run - for aligning raw reads
def initial_bwa_aln(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        bwa_aln(
            args,
            15, 
            ''.join([prefix, '_1.fastq']),
            ''.join([prefix, '_1.sai'])
        )
        bwa_aln(
            args,
            15, 
            ''.join([prefix, '_2.fastq']),
            ''.join([prefix, '_2.sai'])
        )
    else:
        bwa_aln(
            args,
            15, 
            ''.join([prefix, '.fastq']),
            ''.join([prefix, '.sai'])
        )




# Initial SAM creation - for raw read alignment
def initial_bwa_sam(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        bwa_sampe(
            args,
            ''.join([prefix, '_1.fastq']),
            ''.join([prefix, '_1.sai']),
            ''.join([prefix, '_2.fastq']),
            ''.join([prefix, '_2.sai']),
            ''.join([prefix, '.sam'])
        )
        subprocess.call(
            [
                'rm',
                ''.join([prefix, '_1.sai']),
                ''.join([prefix, '_2.sai'])
            ]
        )
    else:
        bwa_samse(
            args,
            ''.join([prefix, '.fastq']),
            ''.join([prefix, '.sai']),
            ''.join([prefix, '.sam'])
        )
    subprocess.call(['rm', ''.join([prefix, '.sai'])])
    



# Initial STAR run - for aligning raw RNA-seq reads
def initial_star(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        star(
            args,
            [
                ''.join([prefix, '_1.fastq']),
                ''.join([prefix, '_2.fastq'])
            ],
            prefix
        )
    else:
        star(
            args,
            [
                ''.join([prefix, '.fastq'])
            ],
            prefix
        )




# Perform the initial sort and index
def initial_sort_and_index(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.bam_start:
        samtools_view(
            args, 
            ''.join([prefix, '.bam']),
            ''.join([prefix, '.filt.bam'])
        )
        samtools_sort(
            args, 
            ''.join([prefix, '.filt.bam']),
            ''.join([prefix, '.sort.bam'])
        )
        subprocess.call(['rm', ''.join([prefix, '.filt.bam'])])
    else:
        if not args.rna_seq:
            samtools_view(
                args, 
                ''.join([prefix, '.sam']),
                ''.join([prefix, '.bam'])
            )
        else:
            samtools_view(
                args, 
                ''.join([prefix, 'Aligned.out.sam']),
                ''.join([prefix, '.bam'])
            )
        samtools_sort(
            args, 
            ''.join([prefix, '.bam']),
            ''.join([prefix, '.sort.bam'])
        )
        subprocess.call(
            [
                'rm',
                ''.join([prefix, '.bam']),
                ''.join([prefix, '.sam'])
            ]
        )
    samtools_index(''.join([prefix, '.sort.bam']))




# Find intersecting SNPs
def find_intersecting_snps(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    subprocess.call(
        [
            args.anaconda_path,
            ''.join(
                [
                    args.wasp_directory,
                    '/mapping/find_intersecting_snps.py'
                ]
            ),
            '--is_sorted',
            '--output_dir', ''.join(
                [
                    args.operating_directory,
                    '/',
                    lane_name
                ]
            ),
            ''.join([prefix, '.sort.bam'])
        ] + [
            '--snp_tab', ''.join([args.hdf5_directory, '/snp_tab.h5']),
            '--snp_index', ''.join([args.hdf5_directory, '/snp_index.h5']),
            '--haplotype', ''.join([args.hdf5_directory, '/haplotypes.h5'])
        ] * (1-args.text_based_snp_files) + [
            '--snp_dir', args.hdf5_directory
        ] * args.text_based_snp_files + [
            '--is_paired_end'
        ] * args.paired_end
    )
    subprocess.call(
        [
            'rm',
            ''.join([prefix, '.sort.bam']),
            ''.join([prefix, '.sort.bam.bai'])
        ]
    )




# Repeat bwa aln for reads overlapping snps
def re_bwa_aln(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        bwa_aln(
            args,
            0, 
            ''.join([prefix, '.sort.remap.fq1.gz']),
            ''.join([prefix, '.remap1.sai'])
        )
        bwa_aln(
            args,
            0, 
            ''.join([prefix, '.sort.remap.fq2.gz']),
            ''.join([prefix, '.remap2.sai'])
        )
    else:
        bwa_aln(
            args,
            0, 
            ''.join([prefix, '.sort.remap.fq.gz']),
            ''.join([prefix, '.remap.sai'])
        )




# Repeat bwa sam for reads overlapping snps
def re_bwa_sam(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        bwa_sampe(
            args,
            ''.join([prefix, '.sort.remap.fq1.gz']),
            ''.join([prefix, '.remap1.sai']),
            ''.join([prefix, '.sort.remap.fq2.gz']),
            ''.join([prefix, '.remap2.sai']),
            ''.join([prefix, '.remap.sam'])
        )
    else:
        bwa_samse(
            args,
            ''.join([prefix, '.sort.remap.fq.gz']),
            ''.join([prefix, '.remap.sai']),
            ''.join([prefix, '.remap.sam'])
        )




# Repeat STAR for reads overlapping snps
def re_star(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        star(
            args,
            [
                ''.join([prefix, '.sort.remap.fq1.gz']),
                ''.join([prefix, '.sort.remap.fq2.gz'])
            ],
            ''.join([prefix, '.remap'])
        )
    else:
        star(
            args,
            [
                ''.join([prefix, '.sort.remap.fq.gz'])
            ],
            ''.join([prefix, '.remap'])
        )




# Remove intermediate files to avoid clutter
def clean_up_fq_sai(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        subprocess.call(
            [
                'rm',
                ''.join([prefix, '.sort.remap.fq1.gz']),
                ''.join([prefix, '.sort.remap.fq2.gz']),
                ''.join([prefix, '.sort.remap.single.fq.gz'])
            ] + [
                ''.join([prefix, '.remap1.sai']),
                ''.join([prefix, '.remap2.sai'])
            ] * (
                1 - args.rna_seq
            )
        )
    else:
        subprocess.call(
            [
                'rm',
                ''.join([prefix, '.sort.remap.fq.gz'])
            ] + [
                ''.join([prefix, '.remap.sai'])
            ] * (
                1 - args.rna_seq
            )
        )




# Convert remapped sam to bam, sort, and index
def re_sort_and_index(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if not args.rna_seq:
        samtools_view(
            args, 
            ''.join([prefix, '.remap.sam']),
            ''.join([prefix, '.remap.bam'])
        )
        subprocess.call(['rm', ''.join([prefix, '.remap.sam'])])
    else:
        samtools_view(
            args, 
            ''.join([prefix, '.remapAligned.out.sam']),
            ''.join([prefix, '.remap.bam'])
        )
        subprocess.call(['rm', ''.join([prefix, '.remapAligned.out.sam'])])
    samtools_sort(
        args, 
        ''.join([prefix, '.remap.bam']),
        ''.join([prefix, '.remap.sort.bam'])
    )
    samtools_index(''.join([prefix, '.remap.sort.bam']))
    subprocess.call(['rm', ''.join([prefix, '.remap.bam'])])




# Filter the remapped reads
def filter_remapped_reads_and_merge(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    subprocess.call(
        [
            args.anaconda_path,
            ''.join([args.wasp_directory, '/mapping/filter_remapped_reads.py']),
            ''.join([prefix, '.sort.to.remap.bam']),
            ''.join([prefix, '.remap.sort.bam']),
            ''.join([prefix, '.remap.keep.bam'])
        ]
    )
    subprocess.call(
        [
            'samtools',
            'merge',
            ''.join([prefix, '.remap.keep.merge.bam']),
            ''.join([prefix, '.sort.keep.bam']),
            ''.join([prefix, '.remap.keep.bam'])
        ]
    )
    subprocess.call(
        [
            'rm',
            ''.join([prefix, '.sort.keep.bam']),
            ''.join([prefix, '.remap.keep.bam']),
            ''.join([prefix, '.sort.to.remap.bam']),
            ''.join([prefix, '.remap.sort.bam']),
            ''.join([prefix, '.remap.sort.bam.bai'])
        ]
    )




# Perform the third sort and index
def sort_and_index_filtered(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    if args.paired_end:
        samtools_sort(args, 
            ''.join([prefix, '.remap.keep.merge.bam']),
            ''.join([prefix, '.remap.keep.merge.sort.sup.bam'])
        )
        samtools_index(''.join([prefix, '.remap.keep.merge.sort.sup.bam']))
        remove_supplementary(args, 
            ''.join([prefix, '.remap.keep.merge.sort.sup.bam']),
            ''.join([prefix, '.remap.keep.merge.sort.bam'])
        )
        subprocess.call(
            [
                'rm',
                ''.join([prefix, '.remap.keep.merge.sort.sup.bam']),
                ''.join([prefix, '.remap.keep.merge.sort.sup.bam.bai'])
            ]
        )
    else:
        samtools_sort(args, 
            ''.join([prefix, '.remap.keep.merge.bam']),
            ''.join([prefix, '.remap.keep.merge.sort.bam'])
        )
    subprocess.call(['rm', ''.join([prefix, '.remap.keep.merge.bam'])])
    samtools_index(''.join([prefix, '.remap.keep.merge.sort.bam']))




# Remove duplicate reads (avoiding bias)
def remove_duplicates(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    subprocess.call(
        [
            args.anaconda_path,
            ''.join(
                [
                    args.wasp_directory
                ] + [
                    '/mapping/rmdup.py'
                ] * (
                    1 - args.paired_end
                ) + [
                    '/mapping/rmdup_pe.py'
                ] * (
                    args.paired_end
                )
            ),
            ''.join([prefix, '.remap.keep.merge.sort.bam']),
            ''.join([prefix, '.remap.keep.merge.rmdup.bam'])
        ]
    )
    subprocess.call(
        [
            'rm',
            ''.join([prefix, '.remap.keep.merge.sort.bam']),
            ''.join([prefix, '.remap.keep.merge.sort.bam.bai'])
        ]
    )




# Perform the final sort and index
def sort_and_index_dedupped(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    samtools_sort(args, 
        ''.join([prefix, '.remap.keep.merge.rmdup.bam']),
        ''.join([prefix, '.remap.keep.merge.rmdup.sort.bam'])
    )
    samtools_index(''.join([prefix, '.remap.keep.merge.rmdup.sort.bam']))
    subprocess.call(['rm', ''.join([prefix, '.remap.keep.merge.rmdup.bam'])])




# Generate a pileup
def pileup(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    subprocess.call(
        [
            'samtools',
            'mpileup',
            '-B',
            '-f', args.reference_genome,
            '-l', args.variant_positions_file,
            '-o', ''.join([prefix, '.pileup']),
            ''.join([prefix, '.remap.keep.merge.rmdup.sort.bam']),
        ]
    )




# Perform binomial tests
def binomial_test(args, lane_name):
    prefix = '/'.join([args.operating_directory, lane_name, lane_name])
    subprocess.call(
        [
            'Rscript',
            '/lab/aaylward-shared/wasp-mapping/binom.r',
            prefix,
            '6',
            str(args.false_discovery_rate),
            args.rsid_directory
        ]
    )
    



#--------------------------- High-level functions -----------------------------#

# Import the relevant lane_names
def create_lane_name_list(args):
    lane_names = []
    with open(args.lane_names) as f:
        for line in f:
            lane_names.append(line.replace('\n', ''))
    return lane_names




# Main
def main(args):

    # Import the list of samples to use
    lane_names = create_lane_name_list(args)
    
    # Construct the arguments that will be passed to pool for multiprocessing
    starmap_arguments = [(args, lane_name) for lane_name in lane_names]
    
    # If the input files are raw reads in fastq format, carry out the initial
    # alignment process
    if not args.bam_start:
        if not args.rna_seq:
            for lane_name in lane_names:
                initial_bwa_aln(args, lane_name)
            pool = Pool(
                processes=min(
                    (
                        args.max_processes, 
                        int(args.memory_limit/5),
                        len(lane_names)
                    )
                )    
            )
            pool.starmap(initial_bwa_sam, starmap_arguments)
            pool.close()
            pool.join()
        else:
            for lane_name in lane_names:
                initial_star(args, lane_name)
    
    # Sort and index the input BAM files, also filtering for quality
    for lane_name in lane_names:
        initial_sort_and_index(args, lane_name)
    
    # Carry out the "find intersecting snps" step
    pool = Pool(
        processes=min(
            (
                args.max_processes, 
                int(args.memory_limit/4),
                len(lane_names)
            )
        )    
    )
    pool.starmap(find_intersecting_snps, starmap_arguments)
    pool.close()
    pool.join()
    
    # Remap indicated reads.
    if not args.rna_seq:
        for lane_name in lane_names:
            re_bwa_aln(args, lane_name)
        pool = Pool(
            processes=min(
                (
                    args.max_processes,
                    int(args.memory_limit/5),
                    len(lane_names)
                )
            )
        )
        pool.starmap(re_bwa_sam, starmap_arguments)
        pool.close()
        pool.join()
    else:
        for lane_name in lane_names:
            re_star(args, lane_name)
        
    # Clean up intermediate files, and sort and index the remapped BAM files.
    for lane_name in lane_names:
        clean_up_fq_sai(args, lane_name)
        re_sort_and_index(args, lane_name)
    
    # Fiter the remapped reads and merge them with the kept reads, sort and
    # index the merged files
    pool = Pool(
        processes=min(
            (
                args.max_processes, 
                int(args.memory_limit),
                len(lane_names)
            )
        )    
    )
    pool.starmap(filter_remapped_reads_and_merge, starmap_arguments)
    pool.close()
    pool.join()
    for lane_name in lane_names:
        sort_and_index_filtered(args, lane_name)
        
    # Remove duplicate reads, sort and index the dedupped files.
    pool = Pool(
        processes=min(
            (
                args.max_processes,
                len(lane_names)
            )
        )    
    )
    pool.starmap(remove_duplicates, starmap_arguments)
    pool.close()
    pool.join()
    for lane_name in lane_names:
        sort_and_index_dedupped(args, lane_name)
    
    # Generate pileups
    pool = Pool(
        processes=min(
            (
                args.max_processes, 
                int(args.memory_limit),
                len(lane_names)
            )
        )    
    )
    pool.starmap(pileup, starmap_arguments)
    pool.close()
    pool.join()
    
    # Perform binomial tests
    pool = Pool(
        processes=min(
            (
                args.max_processes,
                int(args.memory_limit),
                len(lane_names)
            )
        )    
    )
    pool.starmap(binomial_test, starmap_arguments)
    pool.close()
    pool.join()




# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Implementation of the re-mapping procedure to eliminate reference'
            'bias detailed in WASP.'
        )
    )
    parser.add_argument(
        '--anaconda_path',
        required=True,
        type=str,
        help='Path to Anaconda distribution'
    )
    parser.add_argument(
        '--bam_start',
        required=False,
        action='store_true',
        help='Input files are BAM files, skip initial alignment'
    )
    parser.add_argument(
        '--false_discovery_rate',
        required=True,
        type=float,
        help='Set false discovery rate for binomial tests'
    )
    parser.add_argument(
        '--hdf5_directory',
        required=True,
        type=str,
        help='Path to HDF5 snp files directory'
    )
    parser.add_argument(
        '--lane_names',
        required=True,
        type=str,
        help='Path to text file listing lane_names to be considered'
    )
    parser.add_argument(
        '--max_processes',
        required=True,
        type=int,
        help='Maximum number of processes allowed'
    )
    parser.add_argument(
        '--memory_limit',
        required=True,
        type=int,
        help='Approximate memory limit in gigabytes'
    )
    parser.add_argument(
        '--operating_directory',
        required=True,
        type=str,
        help='Path to operating directory'
    )
    parser.add_argument(
        '--paired_end',
        required=False,
        action='store_true',
        help='Input data is paired-end'
    )
    parser.add_argument(
        '--reference_genome',
        required=True,
        type=str,
        help='Path to reference genome'
    )
    parser.add_argument(
        '--rna_seq',
        required=False,
        action='store_true',
        help='Use splice-aware alignment for RNA-seq reads'
    )
    parser.add_argument(
        '--rsid_directory',
        required=True,
        type=str,
        help='Path to RSID directory'
    )
    parser.add_argument(
        '--STAR_genomeDir',
        required=False,
        type=str,
        help='Path to STAR genome directory'
    )
    parser.add_argument(
        '--STAR_sjdbGTFfile',
        required=False,
        type=str,
        help='Path to GTF file'
    )
    parser.add_argument(
        '--text_based_snp_files',
        required=False,
        action='store_true',
        help='use text-based snp files instead'
    )
    parser.add_argument(
        '--variant_positions_file',
        required=True,
        type=str,
        help='Path to variant positions file'
    )
    parser.add_argument(
        '--wasp_directory',
        required=True,
        type=str,
        help='Path to WASP directory'
    )
    return parser.parse_args()




#--------------------------------- Execute ------------------------------------#

args = parse_arguments()
if args.memory_limit < 5:
    raise Exception('Please provide at least 5 GB of memory')
main(args)
