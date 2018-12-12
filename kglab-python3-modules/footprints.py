#!/usr/bin/env python3
#===============================================================================
# footprints.py
#===============================================================================

"""Easy management of transcription factor footprinting data

A mini-module for managing footprint data. The language of this module treats
"footprints" as an abstraction, but mostly handles them as a JSON object / 
dictionary containing the results of a centipede fit.

Example
-------
with footprints.Footprints(
    input_bam='bytes object or path to BAM file',
    input_peaks='bytes object or path to peaks file',
    motifs_db_path='path to motifs database'
) as fp:
    fp.cleans_up_footprints = False
    fp.dump_json('Path to output footprints file')

High-level class
----------------
Footprints
    object representing footprints

Low-level classes
-----------------
Centipede
    class containing Centipede commands
NamedPipe
    context manager for a named pipe

Function
---------
check_input
    check that an input is str or bytes
"""




# Constants ====================================================================

CENTIPEDE_SCRIPT = (
'''#!/usr/bin/R

sink("/dev/null")

library(jsonlite)
library(CENTIPEDE)
library(CENTIPEDE.tutorial)
suppressMessages(library(Rsamtools))

args=commandArgs(trailingOnly=T)
motifs <- args[1]
bam.sorted <- args[2]
test_motif <- args[3]
  
cen <- centipede_data(
    bam_file=bam.sorted,
    fimo_file=motifs,
    pvalue=1e-4,
    flank_size=100
)

fit <- fitCentipede(
  Xlist = list(DNase = cen$mat),
  Y = as.matrix(data.frame(
    Intercept = rep(1, nrow(cen$mat))
  ))
)

json_fit <- toJSON(fit)

sink()

cat(json_fit)
'''
)




# Imports ======================================================================

import json
import os
import os.path
import socket
import subprocess
import sys
import tempfile
from multiprocessing import Pool

hostname = socket.gethostname()
if hostname == 'gatsby.ucsd.edu':
    sys.path.append('/home/data/kglab-python3-modules')
elif hostname == 'holden':
    sys.path.append('/lab/kglab-python3-modules')

import hg19
import namedpipe




# Classes ======================================================================

class Footprints():
    """TF footprints"""
    
    def __init__(
        self,
        input_bam,
        input_peaks,
        motifs_db_path,
        centipede_path=None,
        input_bam_index=None,
        reference_genome_path=hg19.PATH,
        processes=1
    ):
        self.bam = check_input(input_bam)
        self.bam_index = (
            check_input(input_bam_index)
            if
            input_bam_index
            else
            None
        )
        self.peaks = check_input(input_peaks)
        self.motifs_db_path = motifs_db_path
        self.centipede_path = centipede_path
        self.reference_genome_path = reference_genome_path
        self.processes = processes
        self.cleans_up_footprints = False
        self.cleans_up_motifs = False
        self.footprints_file_path = None
        self.motifs_file_path = None
        self.get_fasta()
        self.find_motifs()
        self.infer_footprints()

    def __enter__(self):
        self.cleans_up_footprints = True
        self.cleans_up_motifs = True
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if self.cleans_up_footprints:
            self.clean_up(self.footprints_file_path)
        if self.cleans_up_motifs:
            self.clean_up(self.motifs_file_path)
        return False
    
    def __repr__(self):
        return '\n'.join(
            'Footprints(',
            ')'
        ).format(self)
    
    def get_fasta(self):
        with namedpipe.temp_named_pipe() as peaks_pipe:
            with subprocess.Popen(
                (
                    'sh',
                    '-c',
                    (
                        'bedtools getfasta -bed {0} -fi {1} & cat > {0}'
                        .format(peaks_pipe.name, self.reference_genome_path)
                    )
                ),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE
            ) as bedtools_getfasta:
                self.fasta, _ = bedtools_getfasta.communicate(input=self.peaks)
    
    def find_motifs(self):
        with tempfile.NamedTemporaryFile() as temp_fasta:
            temp_fasta.write(self.fasta)
            with subprocess.Popen(
                (
                    'fimo',
                    '--parse-genomic-coord',
                    '--max-strand',
                    '--text',
                    self.motifs_db_path,
                    temp_fasta.name
                ),
                stdout=subprocess.PIPE
            ) as fimo:
                self.motifs, _ = fimo.communicate()
    
    def infer_footprints(self):
        motif_names = {
            line.split('\t')[0]
            for
            line
            in
            self.motifs.decode().splitlines()
        }
        if self.centipede_path:
            centipede_path = self.centipede_path
        else:
            with tempfile.NamedTemporaryFile() as temp_centipede:
                centipede_path = temp_centipede.name
            with open(centipede_path, 'w') as f:
                f.write(CENTIPEDE_SCRIPT)
        with Pool(processes=min(self.processes, len(motif_names))) as pool:
            self.footprints = dict(
                tup
                for
                tup
                in
                pool.starmap(
                    centipede,
                    (
                        (self, centipede_path, search_motif_name)
                        for
                        search_motif_name
                        in
                        motif_names
                    )
                )
                if
                tup
            )
        os.remove(centipede_path)
    
    def dump_json(self, footprints_file_path):
        with open(footprints_file_path, 'w') as f:
            json.dump(self.footprints, f)
        self.footprints_file_path = footprints_file_path
    
    def write_footprints(self, footprints_file_path):
        self.dump_json(footprints_file_path)
    
    def write_motifs(self, motifs_file_path):
        with open(motifs_file_path, 'wb') as f:
            f.write(self.motifs)
        self.motifs_file_path = motifs_file_path
    
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


def centipede(footprints, centipede_path, search_motif):
    decoded_motifs = footprints.motifs.decode().splitlines()
    header = (
        decoded_motifs[0]
        .replace('sequence_name', 'sequence.name')
        .replace('\tq-value', '')
    )
    with tempfile.NamedTemporaryFile() as (
        temp_motifs_split
    ), tempfile.NamedTemporaryFile() as (
        temp_bam
    ):
        temp_motifs_split.write(
            (
                header
                +
                '\n'
                +
                '\n'.join(
                    line.replace('\t\t', '\t')
                    for
                    line
                    in
                    decoded_motifs[1:]
                    if
                    line.split('\t')[0] == search_motif
                )
                +
                '\n'
            )
            .encode()
        )
        temp_bam.write(footprints.bam)
        if footprints.bam_index:
            with open('{}.bai'.format(temp_bam.name), 'wb') as temp_bai:
                temp_bai.write(footprints.bam_index)
        with subprocess.Popen(
            (
                'Rscript',
                centipede_path,
                temp_motifs_split.name,
                temp_bam.name,
                search_motif
            ),
            stdout=subprocess.PIPE
        ) as centipede:
            fit, _ = centipede.communicate()
        if footprints.bam_index:
            os.remove('{}.bai'.format(temp_bam.name))
    return ((search_motif, json.loads(fit)) if fit else None)
