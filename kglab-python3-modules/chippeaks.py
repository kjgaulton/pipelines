#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#                                chippeaks.py                                   #
#------------------------------------------------------------------------------#

'''Easy management of ChIP-seq peak calling data

A mini-module for managing ChIP-seq peak calling data. The language of this
module treats "ChIP-seq peaks" as an abstraction, but mostly handles them as a
BED file stored in memory.

Example
-------
with chippeaks.ChipPeaks(
    input_bam='bytes object or path to BAM file'
) as cp:
    cp.cleans_up_peaks = False
    cp.blacklist(args.blacklist)
    cp.write('path to output BED-like file')

High-level class
----------------
ChipPeaks
    object representing ChIP-seq peaks

Low-level class
-----------------
NamedPipe
    context manager for a named pipe

Function
---------
check_input
    check that an input is str or bytes
'''




#------------------------------- Dependencies ---------------------------------#

import os
import os.path
import subprocess
import sys
import tempfile

sys.path.append('/home/data/kglab-python3-modules')
sys.path.append('/lab/kglab-python3-modules')

import namedpipe


#---------------------------- Class definitions -------------------------------#

class ChipPeaks():
    '''ChIP-seq peaks'''
    
    def __init__(
        self,
        input_bam
    ):
        self.bam = check_input(input_bam)
        self.cleans_up_peaks = False
        self.peaks_file_path = None
        self.call_peaks()
    
    def __enter__(self):
        self.cleans_up_peaks = True
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.cleans_up_peaks:
            self.clean_up(self.peaks_file_path)
        return False
    
    def __repr__(self):
        return '\n'.join(
            'ChipPeaks(',
            ')'
        ).format(self)
    
    def call_peaks(self):
        with tempfile.NamedTemporaryFile() as (
            temp_bam
        ), tempfile.TemporaryDirectory() as (
            temp_dir_name
        ):
            temp_bam.write(self.bam)
            with tempfile.NamedTemporaryFile(
                dir=temp_dir_name
            ) as temp:
                peaks_pipe_truncated_path = temp.name
            with namedpipe.NamedPipe(
                '{}_peaks.narrowPeak'.format(peaks_pipe_truncated_path)
            ) as peaks_pipe:
                with subprocess.Popen(
                    (
                        'sh',
                        '-c',
                        (
                            'cat {0}_peaks.narrowPeak & '
                            'macs2 callpeak -t {1} -n {0} -B --keep-dup all '
                            '-g hs --extsize 200 --nomodel --shift -100'
                        )
                        .format(
                            peaks_pipe_truncated_path,
                            temp_bam.name
                        )
                    ),
                    stdout=subprocess.PIPE
                ) as macs2:
                    self.peaks, _ = macs2.communicate()
                
    
    def blacklist(self, blacklist_path):
        with subprocess.Popen(
            (
                'bedtools',
                'intersect',
                '-v',
                '-a', 'stdin',
                '-b', blacklist_path,
                '-wa'
            ),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE
        ) as bedtools_intersect:
            peaks, _ = bedtools_intersect.communicate(input=self.peaks)
        self.peaks = peaks
    
    def write(self, peaks_file_path):
        with open(peaks_file_path, 'wb') as f:
            f.write(self.peaks)
        self.peaks_file_path = peaks_file_path
    
    def clean_up(self, path):
        if (os.path.isfile(path) if path else False):
            os.remove(path)




#--------------------------------- Exceptions ---------------------------------#

class Error(Exception):
   '''Base class for other exceptions'''
   pass




class BadInputError(Error):
    '''Bad input error'''
    pass




#---------------------------- Function definitions ----------------------------#

def check_input(input_file):
    '''Check that an input is str or bytes'''
    
    if isinstance(input_file, bytes):
        bytes_obj = input_file
    elif isinstance(input_file, str):
        with open(input_file, 'rb') as f:
            bytes_obj = f.read()
    else:
        raise BadInputError('Input must be either str or bytes')
    return bytes_obj
