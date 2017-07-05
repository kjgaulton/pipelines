#!/usr/bin/env python3

import gzip
import argparse 
import os
import logging
import Levenshtein

def load_adapters(args):
	adapter_files = [os.path.join(args.adapters, f) for f in ['r7_ATAC', 'i7_ATAC', 'i5_ATAC', 'r5_ATAC']]
	for adapter_file in adapter_files:
		if not os.path.exists(adapter_file):
			raise Exception('Missing [{}]. Check if Illumina adapter files [r7_ATAC, i7_ATAC, i5_ATAC, r5_ATAC] exist or are properly named!'.format(os.path.basename(adapter_file)))
	r7_adapters = open(adapter_files[0]).read().splitlines()
	i7_adapters = open(adapter_files[1]).read().splitlines()
	i5_adapters = open(adapter_files[2]).read().splitlines()
	r5_adapters = open(adapter_files[3]).read().splitlines()
	return r7_adapters, i7_adapters, i5_adapters, r5_adapters

def process_mismatch(r7_mismatch, r7_adapters):
	edit_distances = {adap:Levenshtein.distance(r7_mismatch, adap) for adap in r7_adapters}
	min_adap = min(edit_distances, key=edit_distances.get)
	min_dist = edit_distances[min_adap]
	del edit_distances[min_adap]
	min_adap2 = min(edit_distances, key=edit_distances.get)
	min_dist2 = edit_distances[min_adap2]
	return min_adap, min_dist, min_adap2, min_dist2

def correct_barcodes(args, r7_adapters, i7_adapters, i5_adapters, r5_adapters):
	total_reads = 0
	perfect_reads = 0
	skipped_reads = 0
	processed_reads = 0
	
	corrected_output = os.path.join(args.output, args.name + '.correct-bc.fastq.gz')

	with gzip.open(args.reads) as reads_file, gzip.open(corrected_output, 'wb') as corr_out:
		lines = reads_file.read().splitlines()
		for i in range(0, len(lines), 4):
			read_info = ':'.join(lines[i].decode()[1:].split(':')[1:])
			name = lines[i].decode()[1:].split()[0]
			read = lines[i+1].decode()
			plus = lines[i+2].decode()
			qual = lines[i+3].decode()
			barcode = name.split(':')[0]
			r7,i7,i5,r5 = barcode[:8], barcode[8:16], barcode[16:24], barcode[24:]
			total_reads += 1
			if r7 in r7_adapters and i7 in i7_adapters and i5 in i5_adapters and r5 in r5_adapters:
				perfect_reads += 1
			if r7 not in r7_adapters:
				min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(r7, r7_adapters)
				if min_dist <= args.n_mismatches and min_dist2 - min_dist > args.closest_distance:
					r7 = min_adap
				elif min_dist > args.n_mismatches:
					logging.warning('r7 [{}], closest match [{}], edit distance [{}] is greater than maximum n_mismatches [{}]'.format(r7, min_adap, min_dist, args.n_mismatches))
					skipped_reads += 1
					continue
				else:
					logging.warning('r7 [{}], closest match [{}], minimum edit distance [{}] is too close to secondary edit distance [{}]'.format(r7, min_adap, min_dist, min_dist2))
					skipped_reads += 1
					continue

			if i7 not in i7_adapters:
				min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(i7, i7_adapters)
				if min_dist <= args.n_mismatches and min_dist2 - min_dist > args.closest_distance:
					i7 = min_adap
				elif min_dist > args.n_mismatches:
					logging.warning('i7 [{}], closest match [{}], edit distance [{}] is greater than maximum n_mismatches [{}]'.format(i7, min_adap, min_dist, args.n_mismatches))
					skipped_reads += 1
					continue
				else:
					logging.warning('i7 [{}], closest match [{}], minimum edit distance [{}] is too close to secondary edit distance [{}]'.format(i7, min_adap, min_dist, min_dist2))
					skipped_reads += 1
					continue

			if i5 not in i5_adapters:
				min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(i5, i5_adapters)
				if min_dist <= args.n_mismatches and min_dist2 - min_dist > args.closest_distance:
					i5 = min_adap
				elif min_dist > args.n_mismatches:
					logging.warning('i5 [{}], closest match [{}], edit distance [{}] is greater than maximum n_mismatches [{}]'.format(i5, min_adap, min_dist, args.n_mismatches))
					skipped_reads += 1
					continue
				else:
					logging.warning('i5 [{}], closest match [{}], minimum edit distance [{}] is too close to secondary edit distance [{}]'.format(i5, min_adap, min_dist, min_dist2))
					skipped_reads += 1
					continue

			if r5 not in r5_adapters:
				min_adap, min_dist, min_adap2, min_dist2 = process_mismatch(r5, r5_adapters)
				if min_dist <= args.n_mismatches and min_dist2 - min_dist > args.closest_distance:
					r5 = min_adap
				elif min_dist > args.n_mismatches:
					logging.warning('r5 [{}], closest match [{}], edit distance [{}] is greater than maximum n_mismatches [{}]'.format(r5, min_adap, min_dist, args.n_mismatches))
					skipped_reads += 1
					continue
				else:
					logging.warning('r5 [{}], closest match [{}], minimum edit distance [{}] is too close to secondary edit distance [{}]'.format(r5, min_adap, min_dist, min_dist2))
					skipped_reads += 1
					continue
			
			new_barcode = r7 + i7 + i5 + r5
			corr_out.write('@{}:{}\n'.format(new_barcode, read_info).encode())
			corr_out.write('{}\n'.format(read).encode())
			corr_out.write('{}\n'.format(plus).encode())
			corr_out.write('{}\n'.format(qual).encode())
			processed_reads += 1

	return total_reads, perfect_reads, skipped_reads, processed_reads

def main(args):
	logging.info('Starting up.')
	logging.info('Reading in barcoded reads file.')
	
	r7_adapters, i7_adapters, i5_adapters, r5_adapters = load_adapters(args)
	logging.info('Correcting barcoded reads using a maximum edit distance of [{}].'.format(args.n_mismatches))
	total_reads, perfect_reads, skipped_reads, processed_reads = correct_barcodes(args, r7_adapters, i7_adapters, i5_adapters, r5_adapters)
	logging.info('Finishing up.')
	logging.info('Processed barcodes (Ambiguous): {}'.format(total_reads))
	logging.info('Perfect barcodes (Min. distance): {}'.format(perfect_reads))
	logging.info('Skipped barcodes (Min. distance): {}'.format(skipped_reads))
	logging.info('Processed barcodes (Min. distance): {}'.format(processed_reads))
	return

def process_args():
	parser = argparse.ArgumentParser(description='Allow for mismatches for barcodes in a single-cell ATAC-seq experiment.')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-r', '--reads', required=True, type=str, help='Path to barcoded reads fastq file.')
	parser.add_argument('-a', '--adapters', required=False, type=str, default='illumina_adapters/', help='Path to directory that contains Illumina adapter sequences.')
	parser.add_argument('--n_mismatches', required=False, type=int, default=3, help='Maximum number of mismatches.')
	parser.add_argument('--closest_distance', required=False, type=int, default=1, help='Minimum edit distance that secondary match needs to be away from the primary match.')
	return parser.parse_args()


args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
