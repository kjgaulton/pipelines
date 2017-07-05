#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import subprocess
import logging
import traceback
import tempfile
import pandas as pd
from multiprocessing import Pool

def make_windows(args):
	if args.chrom_sizes is not None:
		makewindows_cmd = ['bedtools', 'makewindows', '-g', args.chrom_sizes, '-w', str(200)]
		gzip_cmd = ['gzip', '-c']
		windows_200bp = os.path.join(args.output, 'hg19.200bp_windows.bed.gz')
		with open(windows_200bp, 'w') as f:
			mkwindows = subprocess.Popen(makewindows_cmd, stdout=subprocess.PIPE)
			subprocess.call(gzip_cmd, stdin=mkwindows.stdout, stdout=f)
	return windows_200bp

def genome_coverage(args, bam):
	barcode = os.path.basename(bam).split(args.name+'.')[1].split('.rmdup.bam')[0] 
	genomecov_cmd = ['bedtools', 'genomecov', '-ibam', bam, '-pc', '-bg']
	bed_out = os.path.join(args.output, '.'.join([args.name, barcode, 'PE_coverage.bed']))
	with open(bed_out, 'w') as f:
		genomecov = subprocess.call(genomecov_cmd, stdout=f)
	return bed_out

def make_matrix(args):
	bed_frame = pd.read_csv(args.bed, sep='\t', header=None)
	bed_frame.columns = ['chrom', 'start', 'end']
	for bar in args.barcodes:
		bar_bed = os.path.join(args.output, '.'.join([args.name, bar, 'PE_coverage.bed']))
		intersect_cmd = ['bedtools', 'intersect', '-a', args.bed, '-b', bar_bed, '-wao']
		cut_cmd = ['cut', '-f', '7']
		tmp = tempfile.TemporaryFile('w')
		intersect = subprocess.Popen(intersect_cmd, stdout=subprocess.PIPE)
		subprocess.call(cut_cmd, stdin=intersect.stdout, stdout=tmp)
		tmp.seek(0)
		print(tmp.read())
		tmp.close()
	return

def main(args):
	logging.info('Starting up.')
	try:
		if not os.path.isdir(args.output):
			try:
				os.mkdir(args.output)
			except FileExistsError:
				pass
		if args.bed is None:
			logging.info('BED file unspecified, defaulting to 200 bp windows.')
			args.bed = make_windows(args)
		logging.info('Proceed to calculate genome coverage of input bams using [{}] threads.'.format(args.threads))
		rmdup_bams = glob.glob(os.path.join(args.input, '*.rmdup.bam'))
		args.barcodes = [os.path.basename(bam).split(args.name+'.')[1].split('.rmdup.bam')[0] for bam in rmdup_bams] 
		args_list = [(args, bam) for bam in rmdup_bams]
		#pool = Pool(processes=(args.threads))
		#cov_files = pool.starmap(genome_coverage, args_list)
		#pool.close()
		#pool.join()
		make_matrix(args)
	except Exception as e:
		logging.error(e)
		traceback.print_exc(file=sys.stderr)
		sys.exit(1)
	logging.info('Finishing up.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Create a matrix for read coverage for each sample based on 200 bp windows.')
	parser.add_argument('-i', '--input', required=True, type=str, help='Path to input directory containing rmdup bams')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-b', '--bed', required=False, type=str, help='BED file that will be used to define regions [200bp windows]')
	parser.add_argument('-chrom', '--chrom_sizes', required=False, type=str, default='../misc/hg19.chrom.sizes', help='Path to chromosome size files [../misc/hg19.chrom.sizes]')
	parser.add_argument('-t', '--threads', required=False, type=int, default=16, help='Number of threads for calculating genome coverage')


	return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
