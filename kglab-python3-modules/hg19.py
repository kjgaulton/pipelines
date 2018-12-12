#!/usr/bin/env python3
#===============================================================================
# hg19.py
#===============================================================================

"""Get the coordinates of a variant from its RSID, or an RSID from its
coordinates. 

Examples
--------
rs10_coord = coord('rs10')
print(
    'rs10 is on chromosome {0.chr} at position {0.pos}'
    .format(rs10_coord)
)

rs10_coord_tuple = coord_tuple('rs10')
print(
    'rs10 is on chromosome {} at position {}'
    .format(rs10_coord_tuple[0], rs10_coord_tuple[1])
)

rs_something = rsid(chr=1, pos=10019)
print(
    'The RSID of the variant on chromosome 1 at position 10019 is {}.'
    .format(rs_something)
)

Notes
-----
This module might be a good place to put more utilities later.

coord() returns a fancy object, which is useful for writing readable code. 

coord_tuple() returns a tuple, which is more lightweight and useful for going
fast.

rsid() returns an RSID.

Classes
-------
DataDirectory
    storage class for data directory configuration
Coordinates
    The coordinates of a variant

Functions
---------
coord
    get the coordinates and return them as an object
coord_tuple
    get the coordinates and return them as a tuple
rsid
    get the rsid and return it as a string

Global
------
path
    absolute path to the hg19 reference genome
"""




# Imports ======================================================================

import gzip
import subprocess
import os.path
import socket




# Constants ====================================================================

HOSTNAME = socket.gethostname()

PATH = (
    '/data2/broad-resource-bundle-hg19/ucsc.hg19.fasta'
    if
    HOSTNAME == 'holden'
    else
    '/home/data/broad-resource-bundle-hg19/ucsc.hg19.fasta'
)

SORTED_BY_RSID_FORMAT = (
    '/data2/dbSNP/sorted-by-rsid/{}.bed.gz'
    if
    HOSTNAME == 'holden'
    else
    '/home/data/dbSNP/sorted-by-rsid/{}.bed.gz'
)

SORTED_BY_COORD_PATH = (
    '/data2/dbSNP/dbSNP150.rsid.bed.gz'
    if
    HOSTNAME == 'holden'
    else
    '/home/data/dbSNP/dbSNP150.rsid.bed.gz'
)




# Classes ======================================================================

class Coordinates():
    """The coordinates of a variant"""
    
    def __init__(self, chr, pos):
        self.chr = chr
        self.pos = pos
        self.tuple = chr, pos
    
    def __repr__(self):
        return 'Coordinates(chr={}, pos={})'.format(self.chr, self.pos)


class Variant():
    """The id and coordinates of a variant"""
    
    def __init__(self, id, chr, pos):
        self.id = id
        self.chr = chr
        self.pos = pos
        self.tuple = id, chr, pos
    
    def __repr__(self):
        return 'Variant(id={}, chr={}, pos={})'.format(
            self.id,
            self.chr,
            self.pos
        )




# Functions ====================================================================

def coord(rsid):
    """Get the coordinates and return them as an object"""
    
    chr, pos = coord_tuple(rsid)
    return Coordinates(chr, pos)


def coord_tuple(rsid):
    """Get the coordinates and return them as a tuple"""
    
    with subprocess.Popen(
        ('zcat', SORTED_BY_RSID_FORMAT.format(rsid[:4])),
        stdout=subprocess.PIPE
    ) as zcat:
        with subprocess.Popen(
            (
                'awk',
                '$4=="{}" {{print; exit}}'.format(rsid)
            ),
            stdin=zcat.stdout,
            stdout=subprocess.PIPE
        ) as awk:
            dbsnp_line, _ = awk.communicate()
    try:
        chr, _, pos, _, _, _ = dbsnp_line.decode().split('\t')
    except ValueError:
        raise ValueError(
            '{} was not found in the database'.format(rsid)
        )
    return chr[3:], int(pos)


def rsid(chr, pos):
    """Get the rsid and return it as a string"""
    
    with subprocess.Popen(
        (
            'tabix',
            SORTED_BY_COORD_PATH,
            'chr{0}:{1}-{1}'.format(str(chr).replace('chr', ''), pos)
        ),
        stdout=subprocess.PIPE
    ) as tabix:
        dbsnp_line, _ = tabix.communicate()
    try:
        _, _, _, rsid, *rest = dbsnp_line.decode().split('\t')
    except ValueError:
        raise ValueError(
            'A variant at chromosome {}, position {} was not found in the '
            'database'
            .format(chr, pos)
        )
    return rsid


def range(chr, start, end):
    """Generate all variants within a given genomic range"""
    
    with subprocess.Popen(
        (
            'tabix',
            SORTED_BY_COORD_PATH,
            'chr{0}:{1}-{2}'.format(str(chr).replace('chr', ''), start, end)
        ),
        stdout=subprocess.PIPE
    ) as tabix:
        dbsnp_lines, _ = tabix.communicate()
    for dbsnp_line in dbsnp_lines.decode().splitlines():
        chr, _, pos, rsid, *rest = dbsnp_line.split('\t')
        yield Variant(rsid, chr.replace('chr', ''), int(pos))


def generate_coord_rsid_pairs(file):
    for line in file:
        chr, _, pos, rsid, *alleles = line.split()
        yield (chr.replace('chr', ''), int(pos)), rsid


def coord_rsid_dict():
    """A dictionary containing coord: rsid pairs"""
    
    with gzip.open(SORTED_BY_COORD_PATH, 'rt') as f:
        return dict(generate_coord_rsid_pairs(f))




# test =========================================================================

if __name__ == '__main__':
    rs10_coord = coord('rs10')
    print(
        'rs10 is on chromosome {0.chr} at position {0.pos}'
        .format(rs10_coord)
    )

    rs10_coord_tuple = coord_tuple('rs10')
    print(
        'rs10 is on chromosome {} at position {}'
        .format(rs10_coord_tuple[0], rs10_coord_tuple[1])
    )
    
    rs_something = rsid(chr=1, pos=10019)
    print(
        'The RSID of the variant on chromosome 1 at position 10019 is {}.'
        .format(rs_something)
    )
    
    try:
        coord('rs10a')
    except ValueError:
        print('error was handled')
