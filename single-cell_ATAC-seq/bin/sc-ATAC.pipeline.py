#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import traceback
import logging


def extract_barcodes(args):
	extract_dir = os.path.join(args.output, 'extract-bc')
	if not os.path.isdir(extract_dir):
		try:
			os.mkdir(extract_dir)
		except FileExistsError:
			pass
	extract_p1_cmd = [
			'python3', 'sc-ATAC.extract-bc.py', 
			'-i1', args.index1, '-i2', args.index2, 
			'-r', args.paired1, '-N', str(args.ambiguous_bases),
			'-o', extract_dir, '-n', args.name + '.R1']
	extract_p2_cmd = [
			'python3', 'sc-ATAC.extract-bc.py',
			'-i1', args.index1, '-i2', args.index2, 
			'-r', args.paired2, '-N', str(args.ambiguous_bases),
			'-o', extract_dir, '-n', args.name + '.R2']
	extract_p1_log = os.path.join(extract_dir, '.'.join([args.name, 'R1', 'extract-bc.log']))
	extract_p2_log = os.path.join(extract_dir, '.'.join([args.name, 'R2', 'extract-bc.log']))
	with open(extract_p1_log, 'w') as log1, open(extract_p2_log, 'w') as log2:
		logging.info('--Extracting barcodes for reads 1 file [{}]'.format(args.paired1))
		subprocess.call(extract_p1_cmd, stderr=log1)
		logging.info('--Extracting barcodes for reads 2 file [{}]'.format(args.paired2))
		subprocess.call(extract_p2_cmd, stderr=log2)
	extract_p1 = os.path.join(extract_dir, args.name + '.R1.extract-bc.fastq.gz')
	extract_p2 = os.path.join(extract_dir, args.name + '.R2.extract-bc.fastq.gz')
	return extract_p1, extract_p2

def correct_barcodes(args):
	correct_dir = os.path.join(args.output, 'correct-bc')
	if not os.path.isdir(correct_dir):
		try:
			os.mkdir(correct_dir)
		except FileExistsError:
			pass
	correct_p1_cmd = [
			'python3', 'sc-ATAC.correct-bc.py',
			'-o', correct_dir, '-n', args.name + '.R1',
			'-r', args.extract1, '-a', args.adapters,
			'--n_mismatches', str(args.mismatches),
			'--closest_distance', str(args.closest_distance)]
	correct_p2_cmd = [
			'python3', 'sc-ATAC.correct-bc.py',
			'-o', correct_dir, '-n', args.name + '.R2',
			'-r', args.extract2, '-a', args.adapters,
			'--n_mismatches', str(args.mismatches),
			'--closest_distance', str(args.closest_distance)]
	correct_p1_log = os.path.join(correct_dir, '.'.join([args.name, 'R1', 'correct-bc.log']))
	correct_p2_log = os.path.join(correct_dir, '.'.join([args.name, 'R2', 'correct-bc.log']))
	with open(correct_p1_log, 'w') as log1, open(correct_p2_log, 'w') as log2:
		logging.info('--Correcting barcodes for reads 1 file [{}]'.format(args.extract1))
		subprocess.call(correct_p1_cmd, stderr=log1)
		logging.info('--Correcting barcodes for reads 2 file [{}]'.format(args.extract2))
		subprocess.call(correct_p2_cmd, stderr=log2)
	correct1 = os.path.join(correct_dir, args.name + '.R1.correct-bc.fastq.gz')
	correct2 = os.path.join(correct_dir, args.name + '.R2.correct-bc.fastq.gz')
	return correct1, correct2

def align_reads(args):
	align_dir = os.path.join(args.output, 'aligned_bams')	
	if not os.path.isdir(align_dir):
		try:
			os.mkdir(align_dir)
		except FileExistsError:
			pass
	align_cmd = [
			'python3' , 'sc-ATAC.align_reads.py',
			'-o', align_dir, '-n', args.name,
			'-p1', args.correct1, '-p2', args.correct2,
			'-t', str(args.threads), '-m', str(args.memory),
			'-mapq', str(args.map_quality), '-ref', args.reference
			]
	align_log = os.path.join(align_dir, args.name + '.align.log')
	with open(align_log, 'w') as log:
		subprocess.call(align_cmd, stderr=log)
	compiled_bam = os.path.join(align_dir, args.name + '.compiled.bam') 
	return compiled_bam

def split_cells(args):
	split_dir = os.path.join(args.output, 'split_cells')
	if not os.path.isdir(split_dir):
		try:
			os.mkdir(split_dir)
		except FileExistsError:
			pass
	split_cmd = [
			'python3', 'sc-ATAC.split_cells.py',
			'-i', args.compiled_bam, '-o', split_dir,
			'-n', args.name, '-m', str(args.memory),
			'-min', str(args.min_reads), '-mark_dup', args.picard_mark_dup
			]
	split_log = os.path.join(split_dir, args.name + '.split+rmdup.log')
	with open(split_log, 'w') as log:
		subprocess.call(split_cmd, stderr=log)
	return

def main(args):
	logging.info('Starting up.')
	if not os.path.isdir(args.output):
		try:
			os.mkdir(args.output)
		except FileExistsError:
			pass
		except Exception:
			pass
	try:
		logging.info('Extracting barcodes from index fastq files and inserting them into read names.')
		args.extract1, args.extract2 = extract_barcodes(args)
		if os.path.exists(args.extract1) and os.path.exists(args.extract2):
			logging.info('Correcting barcodes using a maximum of [{}] mismatches and minimum distance of [{}].'.format(args.mismatches, args.closest_distance))
			args.correct1, args.correct2 = correct_barcodes(args)
		else:
			logging.error('Check if [{}] and [{}] exist.'.format(args.extract1, args.extract2))
		if os.path.exists(args.correct1) and os.path.exists(args.correct2):
			logging.info('Aligning reads with BWA-MEM using [{}] threads, [{}]G memory and removing duplicates with Picard MarkDuplicates.'.format(args.threads, args.memory))
			args.compiled_bam = align_reads(args)
		else:
			logging.error('Check if [{}] and [{}] exist.'.format(args.correct1, args.correct2))
		if os.path.exists(args.compiled_bam):
			split_cells(args)
		else:
			logging.error('Check if [{}] exists.'.format(args.compiled_bam))
	except Exception as e:
		logging.error(e)
		traceback.print_exc(file=sys.stderr)
		sys.exit(1)
	logging.info('Finishing up.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Pipeline for single cell ATAC-seq analysis')
	parser.add_argument('-i1', '--index1', required=True, type=str, help='Path to index1 fastq file')
	parser.add_argument('-i2', '--index2', required=True, type=str, help='Path to index2 fastq file')
	parser.add_argument('-p1', '--paired1', required=True, type=str, help='Path to paired-end reads 1 fastq file')
	parser.add_argument('-p2', '--paired2', required=True, type=str, help='Path to paired-end reads 2 fastq file')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-a', '--adapters', required=False, type=str, default='../illumina_adapters/', help='Path to directory that contains Illumina adapters')
	parser.add_argument('-t', '--threads', required=False, type=int, default=8, help='Maximum amount of threads to be used [8]')
	parser.add_argument('-mem', '--memory', required=False, type=int, default=4, help='Maximum amount of memory to be used per thread [4]')
	parser.add_argument('-ref', '--reference', required=True, type=str, help='Reference genome for aligning reads')
	parser.add_argument('-mark_dup', '--picard_mark_dup', required=False, type=str, default='../picard/MarkDuplicates.jar', help='Path to picard\'s MarkDuplicates.jar tool')
	parser.add_argument('--ambiguous_bases', required=False, type=int, default=12, help='Maximum number of ambiguous bases N per barcode [12].')
	parser.add_argument('--mismatches', required=False, type=int, default=3, help='Maximum number of mismatches [3].')
	parser.add_argument('--closest_distance', required=False, type=int, default=1, help='Minimum edit distance that secondary match needs to be away from the primary match [1].')
	parser.add_argument('--min_reads', required=False, type=int, default=500, help='Minimum number of reads for output barcodes [3]')
	parser.add_argument('--map_quality', required=False, type=int, default=10, help='Mapping quality score filter for samtools [10]')
	return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
