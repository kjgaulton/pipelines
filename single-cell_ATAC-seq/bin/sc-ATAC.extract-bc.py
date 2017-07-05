#!/usr/bin/env python3

import gzip
import bz2
import argparse
import os
import logging

def process_barcode(args, i1_file, i2_file, reads_file):
	all_files = {'i1_file': i1_file, 'i2_file': i2_file, 'reads_file': reads_file}
	output_file = os.path.join(args.output, args.name + '.extract-bc.fastq.gz')
	
	total_reads = 0
	skipped_reads = 0
	processed_reads = 0

	with gzip.open(output_file, 'wb') as out: 
		while True:
			full, names, reads, pluses, qualities = [{} for i in range(5)]
			for k,f in all_files.items():
				full[k] = f.readline().strip()[1:].decode()
				names[k] = full[k].split(' ')[0]
				reads[k] = f.readline().strip().decode()
				pluses[k] = f.readline().strip().decode()
				qualities[k] = f.readline().strip().decode()
			if full[k] == reads[k] == pluses[k] == qualities[k] == '':
				# EOF
				return total_reads, skipped_reads, processed_reads
			if '' in names.values():
				raise ValueError('You have a blank read name!')
			if len(list(set(names.values()))) > 1:
				raise ValueError('Your read names are not the same!')
			
			r7 = reads['i1_file'][:8]
			i7 = reads['i1_file'][-8:]
			r5 = reads['i2_file'][:8]
			i5 = reads['i2_file'][-8:]
			barcode = r7 + i7 + r5 + i5
			total_reads += 1

			if barcode.count('N') >= args.N_ambiguous_bases:
				logging.warning('Barcode \"{}\" from read \"{}\" has more than {} N\'s, skipping...'.format(barcode, names['reads_file'], args.N_ambiguous_bases))
				skipped_reads += 1
				continue
			
			out.write('@{}:{}\n'.format(barcode, full['reads_file']).encode())
			out.write('{}\n'.format(reads['reads_file']).encode())
			out.write('+\n'.encode())
			out.write('{}\n'.format(qualities['reads_file']).encode())
			processed_reads += 1
	return total_reads, skipped_reads, processed_reads

def main(args):
	logging.info('Starting up.')
	logging.info('Reading in index and read files.')
	
	read_files = [args.index1, args.index2, args.reads]
	if all(f.endswith('.gz') for f in read_files):
		with gzip.open(args.index1) as i1_file, gzip.open(args.index2) as i2_file, gzip.open(args.reads) as reads_file:
			total_reads, skipped_reads, processed_reads = process_barcode(args, i1_file, i2_file, reads_file)
	elif all(f.endswith('.bz2') for f in read_files):
		with bz2.BZ2File(args.index1) as i1_file, bz2.BZ2File(args.index2) as i2_file, bz2.BZ2File(args.reads) as reads_file:
			total_reads, skipped_reads, processed_reads = process_barcode(args, i1_file, i2_file, reads_file)
	elif all(f.endswith('.fastq') or f.endswith('.fq') for f in read_files):
		with open(args.index1) as i1_file, open(args.index2) as i2_file, open(args.reads) as reads_file:
			total_reads, skipped_reads, processed_reads = process_barcode(args, i1_file, i2_file, reads_file)
	else:
		raise OSError('All file extensions must be the same. Check if all of your input files end with either: .gz, .bz2, .fastq, or .fq.')

	logging.info('Finishing up.')
	logging.info('Total raw reads: {}'.format(total_reads))
	logging.info('Skipped barcodes (Ambiguous): {}'.format(skipped_reads))
	logging.info('Processed barcodes (Ambiguous): {}'.format(processed_reads))
	return

def process_args():
	parser = argparse.ArgumentParser(description='Extract barcodes from a single-cell ATAC-seq experiment.')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-i1', '--index1', required=True, type=str, help='Path to index 1 fastq file.')
	parser.add_argument('-i2', '--index2', required=True, type=str, help='Path to index 2 fastq file.')
	parser.add_argument('-r', '--reads', required=True, type=str, help='Path to reads fastq file.')
	parser.add_argument('-N', '--N_ambiguous_bases', required=False, type=int, default=12, help='Maximum number of ambiguous bases (N).')
	return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
