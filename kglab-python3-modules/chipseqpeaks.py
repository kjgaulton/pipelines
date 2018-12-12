#!/usr/bin/env python3
#===============================================================================
# chipseqpeaks.py
#===============================================================================

"""Easy management of ChIP-seq peak calling data

A mini-module for managing ChIP-seq peak calling data. The language of this
module treats "ChIP-seq peaks" as an abstraction, but mostly handles them as 
MACS2 output stored in memory.

Example
-------
with chipseqpeaks.ChipSeqPeaks(
    input_bam='bytes object or path to BAM file'
) as cp:
    cp.cleans_up = False
    cp.blacklist('path to blacklist')
    cp.write('output prefix')

Class
----------------
ChipSeqPeaks
    object representing ChIP-seq peaks

Function
---------
check_input
    check that an input is str or bytes
"""




# Imports ======================================================================

import os
import os.path
import subprocess
import socket
import sys
import tempfile




# Classes ======================================================================

class ChIPSeqPeaks():
    """ChIP-seq peaks"""
    
    def __init__(
        self,
        treatment_bam,
        control_bam=None,
        qvalue=0.05,
        nomodel=False,
        shift=0,
        broad=False,
        broad_cutoff=0.1,
        log=None
    ):
        self.treatment_bam = check_input(treatment_bam)
        self.control_bam = check_input(control_bam) if control_bam else None
        self.qvalue = qvalue
        self.nomodel = nomodel
        self.shift = shift
        self.broad = broad
        self.broad_cutoff = broad_cutoff
        self.cleans_up = False
        self.cleanup_prefix = None
        self.log = log
        self.output_extensions = (
            (
                'peaks.xls',
                'peaks.narrowPeak',
                'summits.bed',
                'treat_pileup.bdg'
            )
            +
            ('control_lambda.bdg',) * bool(control_bam)
            +
            ('peaks.broadPeak', 'peaks.gappedPeak') * broad
            
        )
        self.call_peaks()
    
    def __enter__(self):
        self.cleans_up = True
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.cleans_up:
            for ext in self.output_extensions:
                self.clean_up('{}_{}'.format(self.cleanup_prefix, ext))
        return False
    
    def __repr__(self):
        return '\n'.join(
            'ChipPeaks(',
            ')'
        ).format(self)
    
    def call_peaks(self):
        with tempfile.NamedTemporaryFile() as (
            temp_treatment_bam
        ), tempfile.NamedTemporaryFile() as (
            temp_control_bam
        ), tempfile.TemporaryDirectory() as (
            temp_dir_name
        ):
            temp_treatment_bam.write(self.treatment_bam)
            if self.control_bam:
                temp_control_bam.write(self.control_bam)
            with tempfile.NamedTemporaryFile(
                dir=temp_dir_name
            ) as temp:
                temp_name = temp.name
            subprocess.call(
                (
                    'macs2', 'callpeak',
                    '-B',
                    '--extsize', '200',
                    '--keep-dup', 'all',
                    '--treatment', temp_treatment_bam.name,
                    '--name', temp_name,
                    '--qvalue', str(self.qvalue),
                    '--shift', str(self.shift),
                )
                +
                (
                    ('--control', temp_control_bam.name)
                    *
                    bool(self.control_bam)
                )
                +
                ('--nomodel',) * self.nomodel
                +
                (
                    ('--broad', '--broad-cutoff', str(self.broad_cutoff))
                    *
                    self.broad
                ),
                stderr=self.log
            )
            for ext in self.output_extensions:
                with subprocess.Popen(
                    ('cat', '{}_{}'.format(temp_name, ext)),
                    stdout=subprocess.PIPE
                ) as cat:
                    output_file, _ = cat.communicate()
                    setattr(self, ext.replace('.', '_'), output_file)
    
    def bdgcmp(self):
        self.output_extensions = self.output_extensions + ('ppois.bdg',)
        with tempfile.NamedTemporaryFile() as (
            temp_treat_pileup
        ), tempfile.NamedTemporaryFile() as (
            temp_control_lambda
        ), tempfile.TemporaryDirectory() as (
            temp_dir_name
        ):
            temp_treat_pileup.write(self.treat_pileup_bdg)
            temp_control_lambda.write(self.control_lambda_bdg)
            with tempfile.NamedTemporaryFile(
                dir=temp_dir_name
            ) as temp:
                temp_name = temp.name
            subprocess.call(
                (
                    'macs2', 'bdgcmp',
                    '-t', temp_treat_pileup.name,
                    '-c', temp_control_lambda.name,
                    '-m', 'ppois',
                    '--o-prefix', temp_name,
                    '-p', '0.00001'
                ),
                stderr=self.log
            )
            with subprocess.Popen(
                ('cat', '{}_ppois.bdg'.format(temp_name)),
                stdout=subprocess.PIPE
            ) as cat:
                self.ppois_bdg, _ = cat.communicate()
    
    def blacklist(self, blacklist_path):
        for peaks in (
            (self.peaks_narrowPeak,)
            +
            (
                (self.peaks_broadPeak, self.peaks_gappedPeak)
                if
                self.broad
                else
                ()
            )
        ): 
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
                stdout=subprocess.PIPE,
                stderr=self.log
            ) as bedtools_intersect:
                peaks, _ = bedtools_intersect.communicate(
                    input=peaks
                )
    
    def write(self, prefix, *extensions):
        for ext in (extensions if extensions else self.output_extensions):
            with open('{}_{}'.format(prefix, ext), 'wb') as f:
                f.write(getattr(self, ext.replace('.', '_')))
        self.cleanup_prefix = prefix
    
    def clean_up(self, path):
        if (os.path.isfile(path) if path else False):
            os.remove(path)




# Exceptions ===================================================================

class Error(Exception):
   """Base class for other exceptions"""
   pass


class BadInputError(Error):
    """Bad input error"""
    pass




# Functions ====================================================================

def check_input(input_file):
    """Check that an input is str or bytes"""
    
    if isinstance(input_file, bytes):
        bytes_obj = input_file
    elif isinstance(input_file, str):
        with open(input_file, 'rb') as f:
            bytes_obj = f.read()
    else:
        raise BadInputError('Input must be either str or bytes')
    return bytes_obj
