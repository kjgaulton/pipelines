#!/usr/bin/env python3

import argparse
import subprocess
import os 
import sys
import reprlib
import logging

#=======================================================#
# trim raw reads with trim_galore!
def trim_galore(args):
	trim_galore_cmd = ['trim_galore', '--fastqc', '-q', '10', '-o', args.output, '--paired', args.paired1, args.paired2]
	with open(os.devnull, 'w') as f:
		subprocess.call(trim_galore_cmd, stderr=f)
	trim_output_1 = os.path.join(args.output, os.path.basename(args.paired1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	trim_output_2 = os.path.join(args.output, os.path.basename(args.paired2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	return trim_output_1, trim_output_2

#=======================================================#
# align reads to the reference genome after throwing 
# out bad quality reads using samtools view.
# remove duplicate reads using samtools rmdup.
def process_reads(args):
	output_prefix = os.path.join(args.output, args.name)
	align_log = output_prefix + '.align.log'
	aligned_bam = output_prefix + '.sort.filt.bam'
	rmdup_bam = output_prefix + '.sort.filt.rmdup.bam'

	bwa_mem_cmd = ['bwa', 'mem', '-t', str(args.processes), args.reference, args.paired1, args.paired2]
	quality_filt_cmd = ['samtools', 'view', '-h', '-F', '1548', '-q', str(args.quality), '-@', str(args.processes), '-']
	mito_filt_cmd = ['grep', '-v', 'chrM']
	samtools_sort_cmd = ['samtools', 'sort', '-m', '{}G'.format(args.memory), '-@', str(args.processes), '-']
	
	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		qual_filt = subprocess.Popen(quality_filt_cmd, stdin=bwa_mem.stdout, stdout=subprocess.PIPE, stderr=log)
		mito_filt = subprocess.Popen(mito_filt_cmd, stdin=qual_filt.stdout, stdout=subprocess.PIPE, stderr=log)
		subprocess.call(samtools_sort_cmd, stdin=mito_filt.stdout, stdout=bam_out, stderr=log)
	
	rmdup_cmd = ['samtools', 'rmdup', aligned_bam, rmdup_bam]
	index_cmd = ['samtools', 'index', rmdup_bam]

	if os.path.exists(aligned_bam) and os.path.getsize(aligned_bam) != 0:
		subprocess.call(rmdup_cmd)
	if os.path.exists(rmdup_bam) and os.path.getsize(rmdup_bam) != 0:
		subprocess.call(index_cmd)
	
	return rmdup_bam

#=======================================================#
# call peaks using macs2
def call_peaks(args, input_bam):
	macs2_log = os.path.join(args.output, '.'.join([args.name, 'macs2_callpeaks.log']))
	macs2_cmd = ['macs2', 'callpeak', 
			'-t', input_bam, 
			'--outdir', args.output, 
			'-n', args.name, 
			'--nomodel', 
			'--shift', '-100', 
			'--extsize', '200', 
			'-B', 
			'--keep-dup', 'all' ]
	with open(macs2_log, 'w') as f:
		subprocess.call(macs2_cmd, stderr=f)
	return

#=======================================================#
# cleanup operations
def cleanup(args):
	read_1_name = os.path.basename(args.paired1).split('.fastq.gz')[0]
	read_2_name = os.path.basename(args.paired2).split('.fastq.gz')[0]
	trim_tmp_1 = os.path.join(args.output, read_1_name + '_trimmed.fq.gz')
	trim_tmp_2 = os.path.join(args.output, read_2_name + '_trimmed.fq.gz')
	fastqc_zip_1 = os.path.join(args.output, read_1_name + '_val_1_fastqc.zip')
	fastqc_zip_2 = os.path.join(args.output, read_2_name + '_val_2_fastqc.zip')
	clean = [trim_tmp_1, trim_tmp_2, fastqc_zip_1, fastqc_zip_2]
	try:
		for f in clean:
			os.remove(f)
	except:
		pass
	return

#=======================================================#
def main(args):
	logging.info('Starting up.')

	if not args.skip_trim:
		try:
			logging.info('Trimming reads with trim_galore.')
			args.paired1, args.paired2 = trim_galore(args)
		except Exception as e:
			print('Check options -p1 and -p2: ' + repr(e), file=sys.stderr)
	if not args.skip_align:
		try:
			logging.info('Aligning reads with bwa and filtering reads with samtools.')
			rmdup_bam = process_reads(args)
		except Exception as e:
			print('Check options -p1 and -p2: ' + repr(e), file=sys.stderr)
	if not args.skip_peaks:
		try:
			logging.info('Calling peaks with MACS2.')
			call_peaks(args, rmdup_bam)
		except Exception as e:
			print(repr(e), file=sys.stderr)
	logging.info('Cleaning up temporary files.')
	cleanup(args)
	logging.info('Finishing up.')
	return

#=======================================================#
# process arguments
def process_args():
	default_ref_genome='/home/joshchiou/references/ucsc.hg19.fasta'

	parser = argparse.ArgumentParser(description='Pipeline for ATAC-seq to trim reads, align them to a reference genome, and then call peaks.')
	parser.add_argument('-p1', '--paired1', required=False, type=str, help='Path to paired reads (1)')
	parser.add_argument('-p2', '--paired2', required=False, type=str, help='Path to paired reads (2)')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory for processed files')
	parser.add_argument('-n', '--name', required=False, type=str, default='sample', help='Output sample name to prepend')
	parser.add_argument('-ref', '--reference', required=False, type=str, default=default_ref_genome, help='Path to reference genome')
	parser.add_argument('-p', '--processes', required=False, type=int, default=4, help='Number of processes to use')
	parser.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum memory per thread for samtools sort')
	parser.add_argument('-q', '--quality', required=False, type=int, default=10, help='Mapping quality cutoff for samtools')
	parser.add_argument('--skip_trim', required=False, action='store_true', default=False, help='Skip adapter trimming step')
	parser.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip read alignment step')
	parser.add_argument('--skip_peaks', required=False, action='store_true', default=False, help='Skip calling peaks step')

	return parser.parse_args()

#=======================================================#

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
