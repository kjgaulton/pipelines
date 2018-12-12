#!/usr/bin/env python3

import sys
import gzip
import argparse
import logging
import Levenshtein

def load_barcodes(args):
	bcs = {}
	with open(args.barcodes) as f:
		for line in f:
			fields = line.rstrip('\n').split('\t')
			index, adapter = fields[0], fields[1]
			bcs[index] = bcs.get(index, []) + [adapter]
	return bcs

def process_mismatch(mismatch, adapters):
	edit_distances = {a:Levenshtein.distance(mismatch, a) for a in adapters}
	min_adap = min(edit_distances, key=edit_distances.get)
	min_dist = edit_distances[min_adap]
	del edit_distances[min_adap]
	min_adap2 = min(edit_distances, key=edit_distances.get)
	min_dist2 = edit_distances[min_adap2]
	return min_adap, min_dist, min_adap2, min_dist2

def demultiplex(args):
	i1 = gzip.open(args.index1, 'rt')
	i2 = gzip.open(args.index2, 'rt')
	r1 = gzip.open(args.reads1, 'rt')
	r2 = gzip.open(args.reads2, 'rt')
	o1 = gzip.open(args.reads1.split('.fastq.gz')[0] + '.demux.fastq.gz', 'wt')
	o2 = gzip.open(args.reads2.split('.fastq.gz')[0] + '.demux.fastq.gz', 'wt')

	edits_made = 0
	while True:
		i1_name = i1.readline().strip()[1:]
		i1_read = i1.readline().strip()
		i1_plus = i1.readline().strip()
		i1_qual = i1.readline().strip()
		
		i2_name = i2.readline().strip()[1:]
		i2_read = i2.readline().strip()
		i2_plus = i2.readline().strip()
		i2_qual = i2.readline().strip()
	
		r1_name = r1.readline().strip()[1:]
		r1_read = r1.readline().strip()
		r1_plus = r1.readline().strip()
		r1_qual = r1.readline().strip()

		r2_name = r2.readline().strip()[1:]
		r2_read = r2.readline().strip()
		r2_plus = r2.readline().strip()
		r2_qual = r2.readline().strip()
		
		if i1_name == '' or i2_name == '' or r1_name == '' or r2_name == '': 
			print('Finished demultiplexing!', file=sys.stderr)
			break
		if not i1_name.split()[0] == i2_name.split()[0] == r1_name.split()[0] == r2_name.split()[0]:
			print('Reads names {} {} {} {} not the same!'.format(i1_name, i2_name, r1_name, r2_name), file = sys.stderr)
			break

		p7 = i1_read[:8]
		i7 = i1_read[-8:]
		i5 = i2_read[:8]
		p5 = i2_read[-8:]		
		if p7 not in args.bcs['p7']:
			min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(p7, args.bcs['p7'])
			if min_dist <= args.mismatch and abs(min_dist2 - min_dist) > args.closest:
				p7 = min_adap
				edits_made += 1
			else:
				continue
		if i7 not in args.bcs['i7']:
			min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(i7, args.bcs['i7'])
			if min_dist <= args.mismatch and abs(min_dist2 - min_dist) > args.closest:
				i7 = min_adap
				edits_made +=1
			else:
				continue
		if i5 not in args.bcs['i5']:
			min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(i5, args.bcs['i5'])
			if min_dist <= args.mismatch and abs(min_dist2 - min_dist) > args.closest:
				i5 = min_adap
				edits_made += 1
			else:
				continue
		if p5 not in args.bcs['p5']:
			min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(p5, args.bcs['p5'])
			if min_dist <= args.mismatch and abs(min_dist2 - min_dist) > args.closest:
				p5 = min_adap
				edits_made += 1
			else:
				continue
		barcode = p7 + i7 + i5 + p5
		o1.write('@' + barcode + ':' + r1_name + '\n')
		o1.write(r1_read + '\n')
		o1.write('+\n')
		o1.write(r1_qual + '\n')
		o2.write('@' + barcode + ':' + r2_name + '\n')
		o2.write(r2_read + '\n')
		o2.write('+\n')
		o2.write(r2_qual + '\n')
	return edits_made

def main(args):
	logging.info('Starting up.')
	args.bcs = load_barcodes(args)
	edits = demultiplex(args)
	print('Total edits made: {}'.format(edits), file=sys.stderr)
	logging.info('Finishing up.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Demultiplex snATAC-seq reads')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-r1', '--reads1', required=True, type=str, help='Path to reads 1 file')
	io_group.add_argument('-r2', '--reads2', required=True, type=str, help='Path to reads 2 file')
	io_group.add_argument('-i1', '--index1', required=True, type=str, help='Path to index 1 file')
	io_group.add_argument('-i2', '--index2', required=True, type=str, help='Path to index 2 file')
	io_group.add_argument('-b', '--barcodes', required=True, type=str, help='Path to barcodes file')
	
	process_group = parser.add_argument_group('Process arguments')
	process_group.add_argument('-m', '--mismatch', required=False, type=int, default=2, help='Maximum number of mismatches per index')
	process_group.add_argument('-c', '--closest', required=False, type=int, default=1, help='Closest edit distance of 2nd best match to the best match')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
