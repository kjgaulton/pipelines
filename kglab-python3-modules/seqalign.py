#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#                                seqalign.py                                   #
#------------------------------------------------------------------------------#

'''Easy management of sequence alignment data

A mini-module for managing sequence alignment data. The language of this module
treats a "sequence alignment" as an abstraction, but mostly handles it as a BAM 
file stored in memory.

Example
-------
with SequenceAlignment(
    input_file_path='path to input bam or fastq file'
) as sa:
    sa.cleans_up_bam = False
    sa.remove_supplementary_alignments()
    sa.samtools_sort(memory_limit=10)
    sa.samtools_index()
    sa.write('path to output BAM file')

Notes
-----
The "input_file" argument should be a string for single-end reads or for
data that is already aligned. For raw paired-end reads, it should be a tuple 
containing two strings giving the paths to the two fasta / fastq files.

High-level class
----------------
SequenceAlignment
    object representing aligned sequencing data

Low-level classes
-----------------
BwaAligner
    commands for running bwa

Functions
---------
file_format_from_extension
    infer the format of a sequencing data file from its extension
median_read_length
    determine the median length of reads in a fasta or fastq file
'''




#------------------------------- Dependencies ---------------------------------#

import gzip
import itertools
import os
import os.path
import subprocess
import math
import socket
import sys
import tempfile
from Bio import SeqIO

hostname = socket.gethostname()
if hostname == 'gatsby.ucsd.edu':
    sys.path.append('/home/data/kglab-python3-modules')
elif hostname == 'holden':
    sys.path.append('/lab/kglab-python3-modules')

import hg19
import namedpipe




#---------------------------- Class definitions -------------------------------#

class SequenceAlignment():
    '''Aligned sequencing data'''
    
    def __init__(
        self,
        input_file,
        phred_quality_score=30,
        processes=1,
        aligner=None,
        log=None
    ):
        self.aligner=aligner
        self.phred_quality_score = int(phred_quality_score)
        self.processes = int(processes)
        self.index = None
        self.bam_file_path = None
        self.cleans_up_bam = False
        self.log = log
        self.parse_input(input_file)
    
    def __enter__(self):
        self.cleans_up_bam = True
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.cleans_up_bam:
            self.clean_up(self.bam_file_path)
            self.clean_up('{}.bai'.format(self.bam_file_path))
        return False
    
    def __repr__(self):
        return '\n'.join(
            'SequenceAlignment(',
            '    phred_quality_score   : {0.phred_quality_score }',
            '    processes             : {0.processes }',
            '    cleans_up_bam         : {0.cleans_up_bam }',
            ')'
        ).format(self)
    
    def parse_input(self, input_file):
        if not isinstance(input_file, (bytes, tuple, str)):
            raise TypeError('input_file must be bytes, tuple, or str')
        elif isinstance(input_file, bytes):
             self.bam = input_file
        elif isinstance(input_file, tuple):
            if len(input_file) != 2:
                raise ValueError(
                    'If input_file_path is a tuple, it must have length 2'
                )
            self.raw_reads_path = input_file
            self.align_reads()
        elif isinstance(input_file, str):
            format = file_format_from_extension(input_file)
            if format in {'fasta', 'fastq'}:
                self.raw_reads_path = input_file
                self.align_reads()
            elif format in {'sam', 'bam'}:
                with subprocess.Popen(
                    (
                        'samtools',
                        'view',
                        '-bhq', str(self.phred_quality_score),
                        '-@', str(self.processes),
                        input_file
                    ),
                    stdout=subprocess.PIPE,
                    stderr=self.log
                ) as samtools_view:
                    self.bam, _ = samtools_view.communicate()
    
    def align_reads(self):
        if not self.aligner:
            self.aligner=BwaAligner()
        self.aligner.align_reads(self)
    
    def remove_supplementary_alignments(self):
        with subprocess.Popen(
            (
				'samtools',
				'view',
				'-bh',
				'-F', '0x800',
				'-@', str(self.processes)
			),
			stdin=subprocess.PIPE,
			stdout=subprocess.PIPE,
			stderr=self.log
        ) as samtools_view:
            bam, _ = samtools_view.communicate(input=self.bam)
        self.bam = bam
    
    def apply_quality_filter(self):
        with subprocess.Popen(
            (
				'samtools',
				'view',
				'-bh',
				'-F', '1548',
				'-q', str(self.phred_quality_score),
				'-@', str(self.processes)
			),
			stdin=subprocess.PIPE,
			stdout=subprocess.PIPE,
			stderr=self.log
        ) as samtools_view:
            bam, _ = samtools_view.communicate(input=self.bam)
        self.bam = bam
    
    def samtools_index(self):
        with namedpipe.temp_named_pipe() as (
            bam_pipe
        ), namedpipe.temp_named_pipe() as (
            index_pipe
        ):
            with subprocess.Popen(
                (
                    'sh',
                    '-c',
                    (
                        'cat {0} & '
                        'samtools index {1} {0} &'
                        'cat > {1}'
                    )
                    .format(
                        index_pipe.name,
                        bam_pipe.name
                    )
                ),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=self.log
            ) as samtools_index:
                self.index, _ = samtools_index.communicate(input=self.bam)
    
    def remove_mitochondrial_reads(self):
        if not self.index:
            raise Exception(
                'use SequenceAlignment.samtools_index() before using '
                'SequenceAlignment.remove_mitochondrial_reads()'
            )
        with tempfile.NamedTemporaryFile() as temp_bam:
            temp_bam.write(self.bam)
            with open('{}.bai'.format(temp_bam.name), 'wb') as f:
                f.write(self.index)
            with subprocess.Popen(
                (
                    'samtools',
                    'view',
                    '-bh',
                    '-@', str(self.processes),
                    temp_bam.name
                )
                +
                tuple(
                    'chr{}'.format(chromosome)
                    for
                    chromosome
                    in
                    (tuple(range(1,23)) + ('X', 'Y'))
                ),
                stdout=subprocess.PIPE,
                stderr=self.log
            ) as samtools_view:
                bam, _ = samtools_view.communicate()
            self.bam = bam
            os.remove('{}.bai'.format(temp_bam.name))
    
    def samtools_sort(self, memory_limit=5):
        if memory_limit < 5:
            raise MemoryLimitError('Please provide at least 5 GB of memory')
        with namedpipe.temp_named_pipe() as sorted_bam_pipe:
            with subprocess.Popen(
                (
                    'sh',
                    '-c',
                    (
                        'cat {0} & '
                        'samtools sort -m {1}M -@ {2} -o {0}'
                    )
                    .format(
                        sorted_bam_pipe.name,
                        int(1024 / self.processes * memory_limit),
                        self.processes
                    )
                ),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=self.log
            ) as samtools_sort:
                bam, _ = samtools_sort.communicate(input=self.bam)
        self.bam = bam
    
    def remove_blacklisted_reads(self, blacklist_path):
        with subprocess.Popen(
            (
                'bedtools',
                'intersect',
                '-abam', 'stdin',
                '-b', blacklist_path,
                '-v'
            ),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self.log
        ) as bedtools_intersect:
            bam, _ = bedtools_intersect.communicate(input=self.bam)
        self.bam = bam
    
    def write(self, bam_file_path):
        with open(bam_file_path, 'wb') as f:
            f.write(self.bam)
        self.bam_file_path = bam_file_path
        if self.index:
            with open('{}.bai'.format(bam_file_path), 'wb') as f:
                f.write(self.index)
    
    def clean_up(self, path):
        if (os.path.isfile(path) if path else False):
            os.remove(path)




class BwaAligner():
    '''A class with methods for calling BWA'''
    
    def __init__(
        self,
        reference_genome_path=hg19.path,
        trim=0,
        max_reads_for_length_check=int(1e6)
    ):
        self.reference_genome_path = reference_genome_path
        self.trims_reads_by = int(trim) if trim else 0
        self.max_reads_for_length_check = max_reads_for_length_check
    
    def __repr__(self):
        return 'BwaAligner()'
    
    def align_reads(self, sequence_alignment):
        self.sa = sequence_alignment
        median_read_length = get_median_read_length(
            self.sa.raw_reads_path,
            self.max_reads_for_length_check
        )
        if median_read_length <= 85:
            self.bwa_aln()
        elif median_read_length > 85:
            self.bwa_mem()
    
    def bwa_aln(self):
        if isinstance(self.sa.raw_reads_path, tuple):
            with namedpipe.temp_named_pipe() as (
                sai_pipe_0
            ), namedpipe.temp_named_pipe() as (
                sai_pipe_1
            ):
                with subprocess.Popen(
                    (
                        'sh',
                        '-c',
                        (
                            'bwa sampe {0} {1} {2} {3} {4} & '
                            'bwa aln -t {5} -q {6} {0} {3} > {1} & '
                            'bwa aln -t {5} -q {6} {0} {4} > {2} & '
                        )
                        .format(
                            self.reference_genome_path,
                            sai_pipe_0.path,
                            sai_pipe_1.path,
                            self.sa.raw_reads_path[0],
                            self.sa.raw_reads_path[1],
                            math.floor(self.sa.processes / 2),
                            self.trims_reads_by
                        )
                    ),
                    stdout=subprocess.PIPE,
                    stderr=self.sa.log
                ) as bwa_aln_sampe:
                    with subprocess.Popen(
                        (
                            'samtools',
                            'view',
                            '-Sbq', str(self.sa.phred_quality_score),
                            '-@', str(self.sa.processes)
                        ),
                        stdin=bwa_aln_sampe.stdout,
                        stdout=subprocess.PIPE,
                        stderr=self.sa.log
                    ) as samtools_view:
                        self.sa.bam, _ = samtools_view.communicate()
        elif isinstance(self.sa.raw_reads_path, str):
            with namedpipe.temp_named_pipe() as sai_pipe:
                with subprocess.Popen(
                    (
                        'sh',
                        '-c',
                        (
                            'bwa samse {0} {1} {2} & '
                            'bwa aln -t {3} -q {4} {0} {2} > {1}; '
                        )
                        .format(
                            self.reference_genome_path,
                            sai_pipe.name,
                            self.sa.raw_reads_path,
                            self.sa.processes,
                            self.trims_reads_by
                        )
                    ),
                    stdout=subprocess.PIPE,
                    stderr=self.sa.log
                ) as bwa_aln_samse:
                    with subprocess.Popen(
                            (
                                'samtools',
                                'view',
                                '-bhq', str(self.sa.phred_quality_score),
                                '-@', str(self.sa.processes)
                            ),
                            stdin=bwa_aln_samse.stdout,
                            stdout=subprocess.PIPE,
                            stderr=self.sa.log
                        ) as samtools_view:
                            self.sa.bam, _ = samtools_view.communicate()
    
    def bwa_mem(self):
        with subprocess.Popen(
            (
                'bwa',
                'mem'
            )
            +
            (
                self.sa.raw_reads_path
                if
                isinstance(self.sa.raw_reads_path, tuple)
                else
                (self.sa.raw_reads_path,)
            ),
            stdout=subprocess.PIPE,
            stderr=self.sa.log
        ) as bwa_mem:
            with subprocess.Popen(
                (
                    'samtools',
                    'view',
                    '-bhq', str(self.sa.phred_quality_score),
                    '-@', str(self.sa.processes)
                ),
                stdin=bwa_mem.stdout,
                stdout=subprocess.PIPE,
                stderr=self.sa.log
            ) as samtools_view:
                self.sa.bam, _ = samtools_view.communicate()




class StarAligner():
    pass




#--------------------------------- Exceptions ---------------------------------#

class Error(Exception):
   '''Base class for other exceptions'''
   pass




class FileExtensionError(Error):
    '''File extension error'''
    pass




class MemoryLimitError(Error):
    '''File extension error'''
    pass




class MissingInputError(Error):
    '''Missing input error'''
    pass




#---------------------------- Function definitions ----------------------------#

def file_format_from_extension(file_path):
    '''Infer the format of a sequencing data file from its extension'''
    
    if (
        (
            file_path.split('.')[-1] in {'fasta', 'fa'}
        )
        or
        (
            file_path.split('.')[-1] == 'gz'
            and
            (
                file_path.split('.')[-2] in {'fasta', 'fa'}
            )
        )
    ):
        format = 'fasta'
    elif (
        (
            file_path.split('.')[-1] in {'fastq', 'fq'}
        )
        or
        (
            file_path.split('.')[-1] == 'gz'
            and
            (
                file_path.split('.')[-2] in {'fastq', 'fq'}
            )
        )
    ):
        format = 'fastq'
    elif file_path.split('.')[-1] in {'sam', 'bam'}:
        format = file_path.split('.')[-1]
    else:
        raise FileExtensionError(
            'Could not parse file extension of {}'
            .format(os.path.basename(file_path))
        )
    return format




def get_median_read_length(raw_reads_paths, number_of_reads):
    '''Return the median read length of a FASTA or FASTQ file'''
    
    histogram = {}
    if isinstance(raw_reads_paths, tuple):
        formats = tuple(
            file_format_from_extension(raw_reads_path[i])
            for
            i
            in
            range(2)
        )
    elif isinstance(raw_reads_paths, str):
        formats = (file_format_from_extension(raw_reads_paths),)
        raw_reads_paths = (raw_reads_paths,)
    for raw_reads_path, format in zip(raw_reads_paths, formats):
        with (
            gzip.open(raw_reads_path, 'rt')
            if
            raw_reads_path[-3:] == '.gz'
            else
            open(raw_reads_path, 'r')
        ) as raw_reads:
            for record in itertools.islice(
                SeqIO.parse(raw_reads, format),
                number_of_reads
            ):
                try:
                    histogram[len(record.seq)] += 1
                except KeyError:
                    histogram[len(record.seq)] = 1
        read_lengths = tuple(
                    length for length, count in sorted(histogram.items())
                )
        total_reads = sum(count for length, count in histogram.items())
        cumulative_count = 0
        for length, count in sorted(histogram.items()):
            cumulative_count += count
            if cumulative_count > total_reads / 2:
                median = length
                break
            elif cumulative_count == total_reads / 2:
                read_lengths = tuple(
                    length for length, count in sorted(histogram.items())
                )
                next_length = read_lengths[read_lengths.index(length) + 1]
                median = (length + next_length) / 2
                break
    return median




#----------------------------------- Test -------------------------------------#

if __name__ == '__main__':
    with SequenceAlignment(
        input_file='/home/data/aaylward-data/bam/ChIP_FOXA1_hepg2_2.fq',
        processes=23,
        aligner=BwaAligner(
            trim=15
        )
    ) as sa:
        sa.cleans_up_bam=False
        print('applying quality filter')
        sa.apply_quality_filter()
        print('sorting')
        sa.samtools_sort(memory_limit=64)
        print('indexing')
        sa.samtools_index()
        print('removing mitochondrial reads')
        sa.remove_mitochondrial_reads()
        sa.write('/home/data/aaylward-data/bam/joshs-chip-test.bam')
