#!/usr/bin/env python3

import argparse
import subprocess
from pyfaidx import Fasta

def process_args():
	parser = argparse.ArgumentParser(description='Script to insert SNP alleles into reference sequences')
	parser.add_argument('-p', '--peaks', required=True, type=str, help='Path to peaks file (3 col)')
	parser.add_argument('-v', '--variants', required=True, type=str, help='Path to variants bed-like file (6 col)')
	parser.add_argument('-f', '--fasta', required=True, type=str, help='Path to reference fasta')
	return parser.parse_args()

def get_alt_sequence(args, fasta, chrom, start, stop, var, ref, alt):
	header = '>{0}:{1}-{2}\t{0}:{3}_{4}/{5}'.format(chrom, start, stop, var, ref, alt)
	left = fasta[chrom][start:var-1]
	right = fasta[chrom][var:stop]
	ref_sequence = left + ref + right
	alt_sequence = left + alt + right
	return header, ref_sequence, alt_sequence

def load_fasta(args):
	fasta = Fasta(args.fasta, as_raw=True)
	return fasta

def get_atac_variants(args):
	intersect_out = '.'.join([args.peaks.split('.bed')[0], args.variants.split('.bed')[0], 'tmp'])
	cmd = ['bedtools', 'intersect',
			'-a', args.peaks,
			'-b', args.variants,
			'-wa','-wb']
	with open(intersect_out, 'w') as f:
		subprocess.call(cmd, stdout=f)
	return intersect_out

def process_intersect(line):
	fields = line.split('\t')
	chrom, start, stop = fields[0], int(fields[1]), int(fields[2])
	var, ref, alt = int(fields[5]), fields[7], fields[8]
	return chrom, start, stop, var, ref, alt

def main(args):
	intersect_out = get_atac_variants(args)
	fasta = load_fasta(args)
	with open(intersect_out, 'r') as f:
		lines = f.read().splitlines()
		for line in lines:
			chrom, start, stop, var, ref, alts = process_intersect(line)
			for alt in alts.split(','):
				header, ref_sequence, alt_sequence = get_alt_sequence(args, fasta, chrom, start, stop, var, ref, alt)
				print(header)
				print(alt_sequence)
	return

args = process_args()
main(args)
