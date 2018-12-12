#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import logging
import subprocess
import pysam
import numpy as np
import pandas as pd
import scipy.sparse

def trim_reads(args):
	trim_galore_cmd = ['trim_galore', '--nextera', '--fastqc', '-q', '20', '-o', args.output, '--paired', args.read1, args.read2]
	with open(os.devnull, 'w') as f:
		subprocess.call(trim_galore_cmd, stderr=f, stdout=f)
	pair1_trim = os.path.join(args.output, os.path.basename(args.read1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	pair2_trim = os.path.join(args.output, os.path.basename(args.read2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	return pair1_trim, pair2_trim

def align_reads(args):
	align_log = args.output_prefix + '.align.log'
	aligned_bam = args.output_prefix + '.compiled.bam'
	filtered_bam = args.output_prefix + '.compiled.filt.bam'

	bwa_mem_cmd = ['bwa', 'mem', '-t', str(args.threads), args.reference, args.read1, args.read2]
	fixmate_cmd = ['samtools', 'fixmate', '-r', '-', '-']
	samtools_sort_cmd = ['samtools', 'sort', '-m', str(args.memory)+'G', '-@', str(args.threads), '-']
	filter_cmd = ['samtools', 'view', '-h', '-q', str(args.map_quality), '-f', '3', '-F', '4', '-F', '256', '-F', '1024', '-F', '2048', aligned_bam]
	samtools_view_cmd = ['samtools', 'view', '-b', '-']
	
	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		fixmate = subprocess.Popen(fixmate_cmd, stdin=bwa_mem.stdout, stdout=subprocess.PIPE)
		subprocess.call(samtools_sort_cmd, stdin=fixmate.stdout, stdout=bam_out, stderr=log)
	
	with open(align_log, 'a') as log, open(filtered_bam, 'w') as bam_out:
		filt = subprocess.Popen(filter_cmd, stderr=log, stdout=subprocess.PIPE)
		view = subprocess.Popen(samtools_view_cmd, stdin=subprocess.PIPE, stdout=bam_out)
		for line in filt.stdout:
			line = line.decode().rstrip('\n')
			if line.startswith('@'):
				try:
					view.stdin.write('{}\n'.format(line).encode())
				except IOError as err:
					if err.errno == errno.EPIPE or err.errno == errno.EINVAL:
						break
					else:
						raise
				continue
			fields = line.split('\t')
			barcode = fields[0].split(':')[0]
			fields[0] = barcode + '_' + ':'.join(fields[0].split(':')[1:])
			fields.append('BX:Z:{}'.format(barcode))
			try:
				view.stdin.write('{}\n'.format('\t'.join(fields)).encode())
			except IOError as err:
				if err.errno == errno.EPIPE or err.errno == errno.EINVAL:
					break
				else:
					raise
		view.stdin.close()
		view.wait()
	return

def remove_duplicate_reads(args):
	filtered_bam = args.output_prefix + '.compiled.filt.bam'
	markdup_bam = args.output_prefix + '.filt.md.bam'
	rmdup_bam = args.output_prefix + '.filt.rmdup.bam'
	markdup_log = args.output_prefix + '.MarkDuplicates.log'
	markdup_cmd = ['java', '-Xmx24G', '-jar', args.picard, 
			'MarkDuplicates', 'INPUT={}'.format(filtered_bam), 'OUTPUT={}'.format(markdup_bam), 
			'VALIDATION_STRINGENCY=LENIENT', 'BARCODE_TAG=BX', 'METRICS_FILE={}'.format(markdup_log),
			'REMOVE_DUPLICATES=false']
	index_cmd = ['samtools', 'index', markdup_bam]
	filter_cmd = ['samtools', 'view', '-@', str(args.threads), '-b', '-f', '3', '-F', '1024', markdup_bam]
	filter_cmd.extend(['chr{}'.format(c) for c in list(map(str, range(1,23))) + ['X','Y']])
	
	with open(os.devnull, 'w') as null:
		subprocess.call(markdup_cmd, stderr=null, stdout=null)
		subprocess.call(index_cmd)	
	if os.path.isfile(markdup_bam):
		with open(rmdup_bam, 'w') as bam_out:
			subprocess.call(filter_cmd, stdout=bam_out)
	else:
		raise FileNotFoundError('{} not found!'.format(markdup_bam))
	return

def generate_binary_matrix(args):
	rmdup_bam = args.output_prefix + '.filt.rmdup.bam'
	tagalign_file = args.output_prefix + '.filt.rmdup.tagAlign.gz'
	
	barcode_coverage = {}
	if os.path.isfile(rmdup_bam):
		with gzip.open(tagalign_file, 'wt') as f:
			temp_store = []
			bamfile = pysam.AlignmentFile(rmdup_bam, 'rb')
			genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
			for read in bamfile:
				if not read.is_proper_pair:
					continue
				barcode = read.query_name.split('_')[0]
				read_chr = read.reference_name
				read_start = max(1, read.reference_start - args.extension if read.is_reverse else read.reference_start + 4 - args.extension)
				read_end = min(genome_size[read_chr], read.reference_end - 5 + args.extension if read.is_reverse else read.reference_end + args.extension)
				read_qual = read.mapping_quality
				if read.is_reverse:
					read_orient = '-'
				else:
					read_orient = '+'
				barcode_coverage[barcode] = barcode_coverage.get(barcode, 0) + 1
				line_out = '\t'.join([read_chr, str(read_start), str(read_end), '{}_{}'.format(args.name,barcode), str(read_qual), read_orient])
				print(line_out, sep='\t', file=f)
			bamfile.close()
	else:
		raise FileNotFoundError('{} not found!'.format(rmdup_bam))
	
	barcode_counts = args.output_prefix + '.barcode_counts.txt'
	with open(barcode_counts, 'w') as f:
		print('barcode', 'count', sep='\t', file=f)
		for bc in barcode_coverage:
			print('{}_{}'.format(args.name,bc), barcode_coverage[bc], sep='\t', file=f)
	pass_barcodes = ['{}_{}'.format(args.name, bc) for bc in barcode_coverage if barcode_coverage[bc] >= args.minimum_reads]

	barcodes_file = args.output_prefix + '.barcodes'
	peaks_file = args.output_prefix + '.peaks'
	mtx_file = args.output_prefix + '.int.csr.npz'

	matrix = {}

	regions_file = make_windows(args)
	peak_intersect = intersect_helper(tagalign_file, regions_file)

	for line in peak_intersect.stdout:
		line = line.decode().rstrip('\n')
		fields = line.split('\t')
		barcode = fields[3]
		peak = fields[9]
		if barcode not in pass_barcodes:
			continue
		if peak not in matrix:
			matrix[peak] = {}
		matrix[peak][barcode] = matrix[peak].get(barcode, 0) + 1

	threshold = 2
	for peak in list(matrix):
		if len(matrix[peak].keys()) < threshold:
			del matrix[peak]
	peaks = list(matrix)
	barcodes = pass_barcodes
	mtx = scipy.sparse.dok_matrix((len(barcodes), len(peaks)), dtype=int)
	for b in range(len(barcodes)):
		for p in range(len(peaks)):
			if barcodes[b] in matrix[peaks[p]]:
				mtx[b,p] = matrix[peaks[p]][barcodes[b]]
	mtx = mtx.tocsr()
	scipy.sparse.save_npz(mtx_file, mtx)
	
	with open(barcodes_file, 'w') as f:
		print('\n'.join(barcodes), file=f)
	with open(peaks_file, 'w') as f:
		print('\n'.join(peaks), file=f)
	return

def intersect_helper(tagalign, regions):
	awk_cmd = ['awk', '''BEGIN{{FS=OFS="\\t"}} {{peakid=$1":"$2"-"$3; gsub("chr","",peakid); print $1, $2, $3, peakid}}''', regions]
	intersect_cmd = ['bedtools', 'intersect', '-a', tagalign, '-b', '-', '-wa', '-wb']
	awk = subprocess.Popen(awk_cmd, stdout=subprocess.PIPE)
	intersect = subprocess.Popen(intersect_cmd, stdin=awk.stdout, stdout=subprocess.PIPE)
	return intersect

def make_windows(args):
	makewindows_cmd = ['bedtools', 'makewindows', '-g', args.chrom_sizes, '-w', str(args.window_size * 1000)]
	filter_cmd = ['grep', '-v', '_']
	blacklist_cmd = ['bedtools', 'intersect', '-a', '-', '-b', args.blacklist, '-v']
	
	regions_file = '{}.{}kb_windows.bed'.format(args.output_prefix, args.window_size)
	with open(regions_file, 'w') as f:
		makewindows = subprocess.Popen(makewindows_cmd, stdout=subprocess.PIPE)
		filt = subprocess.Popen(filter_cmd, stdin=makewindows.stdout, stdout=subprocess.PIPE)
		subprocess.call(blacklist_cmd, stdin=filt.stdout, stdout=f)
	return regions_file

def main(args):
	logging.info('Start.')
	if not os.path.isdir(args.output):
		os.makedirs(args.output)
	args.output_prefix = os.path.join(args.output, args.name)
	if not args.skip_trim:
		logging.info('Trimming adapters using trim_galore.')
		args.read1, args.read2 = trim_reads(args)
	if not args.skip_align:
		logging.info('Aligning reads using BWA mem with [{}] processors.'.format(args.threads))
		logging.info('Piping output to samtools using [{}] Gb of memory.'.format(args.memory))
		align_reads(args)
	if not args.skip_rmdup:
		logging.info('Removing duplicate and mitochrondrial reads.'.format(args.minimum_reads))
		remove_duplicate_reads(args)
	if not args.skip_matrix:
		logging.info('Generating tagalign and chromatin accessibility matrix.')
		generate_binary_matrix(args)
	logging.info('Finish.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Align demulitplexed snATAC-seq reads to a reference genome.')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-r1', '--read1', required=True, type=str, help='Paired-end reads file 1')
	io_group.add_argument('-r2', '--read2', required=True, type=str, help='Paired-end reads file 2')
	io_group.add_argument('-o', '--output', required=False, type=str, default=os.getcwd(), help='Output directory to store processed files')
	io_group.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	
	align_group = parser.add_argument_group('Alignment arguments')
	align_group.add_argument('-t', '--threads', required=False, type=int, default=8, help='Number of threads to use for alignment [8]')
	align_group.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum amount of memory (G) per thread for samtools sort [4]')
	align_group.add_argument('-q', '--map_quality', required=False, type=int, default=30, help='Mapping quality score filter for samtools [30]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='/home/joshchiou/references/male.hg19.fa', help='Path to the BWA-prepared reference genome')
	
	dup_group = parser.add_argument_group('Remove duplicates arguments')
	dup_group.add_argument('--picard', required=False, type=str, default='/home/joshchiou/bin/picard.jar', help='Path to picard.jar')
	
	matrix_group = parser.add_argument_group('Matrix generation arguments')
	matrix_group.add_argument('--extension', required=False, type=int, default=75, help='Read extension length')
	matrix_group.add_argument('--minimum-reads', required=False, type=int, default=1000, help='Minimum number of reads for barcode inclusion')
	matrix_group.add_argument('--window-size', required=False, type=int, default=5, help='Size (kb) to use for defining windows of accessibility')
	matrix_group.add_argument('--chrom-sizes', required=False, type=str, default='/home/joshchiou/references/hg19.chrom.sizes', help='Chromosome sizes file from UCSC')
	matrix_group.add_argument('--blacklist', required=False, type=str, default='/home/joshchiou/references/ENCODE.hg19.blacklist.bed', help='BED file of blacklisted regions')

	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-trim', required=False, action='store_true', default=False, help='Skip adapter trimming step')
	skip_group.add_argument('--skip-align', required=False, action='store_true', default=False, help='Skip read alignment step')
	skip_group.add_argument('--skip-rmdup', required=False, action='store_true', default=False, help='Skip duplicate removal step')
	skip_group.add_argument('--skip-matrix', required=False, action='store_true', default=False, help='Skip matrix generation step')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
	args = process_args()
	main(args)
