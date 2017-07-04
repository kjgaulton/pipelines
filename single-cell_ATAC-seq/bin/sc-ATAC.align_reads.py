#!/usr/bin/env python3

import os
import argparse
import logging
import subprocess


def trim_reads(args):
	trim_galore_cmd = ['trim_galore', '--nextera', '--fastqc', '-q', '20', '-o', args.output, '--paired', args.paired1, args.paired2]
	with open(os.devnull, 'w') as f:
		subprocess.call(trim_galore_cmd, stderr=f, stdout=f)
	pair1_trim = os.path.join(args.output, os.path.basename(args.paired1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	pair2_trim = os.path.join(args.output, os.path.basename(args.paired2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	return pair1_trim, pair2_trim

def align_reads(args):
	output_prefix = os.path.join(args.output, args.name) 

	bwa_mem_cmd = ['bwa', 'mem', '-t', str(args.threads), args.reference, args.paired1, args.paired2]
	quality_filt_cmd = ['samtools', 'view', '-h', '-f', '2', '-q', str(args.map_quality), '-@', '1', '-']
	mito_filt_cmd = ['grep', '-v', 'chrM']
	samtools_sort_cmd = ['samtools', 'sort', '-m', str(args.memory)+'G', '-@', str(args.threads), '-']
	
	align_log = output_prefix + '.align.log'
	aligned_bam = output_prefix + '.compiled.bam'

	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		qual_filt = subprocess.Popen(quality_filt_cmd, stdin=bwa_mem.stdout, stdout=subprocess.PIPE, stderr=log)
		mito_filt = subprocess.Popen(mito_filt_cmd, stdin=qual_filt.stdout, stdout=subprocess.PIPE, stderr=log)
		subprocess.call(samtools_sort_cmd, stdin=mito_filt.stdout, stdout=bam_out, stderr=log)
	
	count_reads_cmd = ['samtools', 'view', aligned_bam]
	p = subprocess.Popen(count_reads_cmd, stdout=subprocess.PIPE)
	aligned_reads = sum(1 for _ in p.stdout)
	return aligned_reads / 2

def main(args):
	logging.info('Starting up.')
	
	if not args.skip_trim:
		logging.info('Trimming adapters using trim_galore.')
		args.paired1, args.paired2 = trim_reads(args)
	
	if not args.skip_align:
		logging.info('Aligning reads using BWA mem with [{}] processors.'.format(args.threads))
		logging.info('Piping output to samtools using [{}] Gb of memory.'.format(args.memory))
		aligned_reads = align_reads(args)
	logging.info('Finishing up.')
	if aligned_reads is not None:
		logging.info('Aligned pairs: {}'.format(aligned_reads))
	return

def process_args():
	parser = argparse.ArgumentParser(description='Align all reads to the reference genome and make compiled bam.')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-p1', '--paired1', required=True, type=str, help='Paired reads file 1')
	parser.add_argument('-p2', '--paired2', required=True, type=str, help='Paired reads file 2')
	parser.add_argument('-t', '--threads', required=False, type=int, default=8, help='Number of threads to use for alignment')
	parser.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum amount of memory to use for samtools sort')
	parser.add_argument('-mapq', '--map_quality', required=False, type=int, default=10, help='Mapping quality score filter for samtools')
	parser.add_argument('-ref', '--reference', required=True, type=str, help='Path to the reference genome')
	parser.add_argument('--skip_trim', required=False, action='store_true', default=False, help='Skip adapter trimming step')
	parser.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip aligning reads step')
	return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
