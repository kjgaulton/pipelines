#!/usr/bin/env python3

import argparse
import subprocess
import os 
import sys
import reprlib
import logging

#=======================================================#

def process_reads(args, reads, name):
	output_prefix = os.path.join(args.output, name)
	align_log = output_prefix + '.align.log'
	aligned_bam = output_prefix + '.sort.filt.bam'
	rmdup_bam = output_prefix + '.sort.filt.rmdup.bam'
	
	# align reads with bwa aln / bwa samse
	# then filter out reads that are unmapped, chrM reads, < quality score
	# then sort reads
	bwa_aln_cmd = ['bwa', 'aln', '-q', '15', '-t', str(args.processes), args.reference, reads]
	bwa_samse_cmd = ['bwa', 'samse', args.reference, '-', reads]
	quality_filt_cmd = ['samtools', 'view', '-h', '-F', '1548', '-q', str(args.quality), '-@', str(args.processes), '-']
	mito_filt_cmd = ['grep', '-v', 'chrM']
	samtools_sort_cmd = ['samtools', 'sort', '-m', '{}G'.format(args.memory), '-@', str(args.processes), '-']
	
	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_aln = subprocess.Popen(bwa_aln_cmd, stdout=subprocess.PIPE, stderr=log)
		bwa_samse = subprocess.Popen(bwa_samse_cmd, stdin=bwa_aln.stdout, stdout=subprocess.PIPE, stderr=log)
		qual_filt = subprocess.Popen(quality_filt_cmd, stdin=bwa_samse.stdout, stdout=subprocess.PIPE, stderr=log)
		mito_filt = subprocess.Popen(mito_filt_cmd, stdin=qual_filt.stdout, stdout=subprocess.PIPE, stderr=log)
		subprocess.call(samtools_sort_cmd, stdin=mito_filt.stdout, stdout=bam_out, stderr=log)
	
	# use picard MarkDuplicates to filter out duplicate reads
	# then index the bam file
	metrics_file = output_prefix + '.MarkDuplicates.metrics'
	rmdup_cmd = [
			'java', '-Xmx{}G'.format(args.memory), 
			'-jar', args.markdup,
			'INPUT={}'.format(aligned_bam),
			'OUTPUT={}'.format(rmdup_bam),
			'REMOVE_DUPLICATES=true',
			'VALIDATION_STRINGENCY=LENIENT',
			'METRICS_FILE={}'.format(metrics_file)]
	index_cmd = ['samtools', 'index', rmdup_bam]

	if os.path.exists(aligned_bam) and os.path.getsize(aligned_bam) != 0:
		with open(os.devnull, 'w') as f:
			subprocess.call(rmdup_cmd, stderr=f)
	if os.path.exists(rmdup_bam) and os.path.getsize(rmdup_bam) != 0:
		subprocess.call(index_cmd)
	
	return rmdup_bam

#=======================================================#

def call_peaks(args, treat_bam, control_bam):
	macs2_log = os.path.join(args.output, '.'.join([args.name, 'macs2_callpeaks.log']))
	macs2_cmd = ['macs2', 'callpeak', 
			'-t', treat_bam,
			'-c', control_bam,
			'--outdir', args.output, 
			'-n', args.name, 
			'--extsize', '200', 
			'-B', 
			'--keep-dup', 'all' ]
	if args.broad:
		macs2_cmd.extend(['--broad'])
	with open(macs2_log, 'w') as f:
		subprocess.call(macs2_cmd, stderr=f)
	return

#=======================================================#

def bdgcmp(args):
	bdgcmp_log = os.path.join(args.output, '.'.join([args.name, 'bdgcmp.log']))
	bdgcmp_cmd = [
			'macs2', 'bdgcmp',
			'-t', os.path.join(args.output, args.name + '_treat_pileup.bdg'),
			'-c', os.path.join(args.output, args.name + '_control_lambda.bdg'),
			'-m', 'ppois',
			'--outdir', args.output,
			'--o-prefix', args.name,
			'-p', '0.00001']
	with open(bdgcmp_log, 'w') as f:
		subprocess.call(bdgcmp_cmd, stderr=f)
	
	bdgcmp_out = os.path.join(args.output, args.name + '_ppois.bdg')
	
	# add track label to the bedgraph, sort, and gzip
	label = 'track type=bedGraph name=\"{0}\" description=\"{0}\" visibility=2 color=0,0,0 altColor=0,0,0 autoScale=off maxHeightPixels=64:64:32'.format(args.name)
	sorted_bdg = os.path.join(args.output, args.name  + '_ppois.sorted.bdg'
	sort_cmd = ['sort', '-k', '1,1', '-k', '2,2n', bdg]
	with open(sorted_bdg, 'w') as f:
		print(label, file=f)
	with open(sorted_bdg, 'a') as f:
		subprocess.call(sort_cmd, stdout=f)
	subprocess.call(['gzip', sorted_bdg])
	os.remove(bdgcmp_out)
	return

#=======================================================#

def main(args):
	logging.info('Starting up.')

	if not args.skip_align:
		try:
			logging.info('Aligning reads with bwa aln/samse and filtering reads with samtools.')
			logging.info('Processing treatment bam [{}].'.format(args.treatment))
			treat_bam = process_reads(args, args.treatment, args.name)
			logging.info('Processing control bam [{}].'.format(args.control))
			control_bam = process_reads(args, args.control, args.name + '_control')
		except Exception as e:
			print('Check options -t and -c: ' + repr(e), file=sys.stderr)
	if not args.skip_peaks:
		try:
			logging.info('Calling peaks with MACS2.')
			call_peaks(args, treat_bam, control_bam)
		except Exception as e:
			print(repr(e), file=sys.stderr)
	if not args.skip_track:
		try:
			logging.info('Building signal track with MACS2 output.')
			bdgcmp(args)
		except Exception as e:
			print(repr(e), file=sys.stderr)
	logging.info('Finishing up.')
	return

#=======================================================#

def process_args():
	parser = argparse.ArgumentParser(description='Pipeline for ChIP to align reads to a reference genome, and then call peaks.')
	parser.add_argument('-t', '--treatment', required=True, type=str, help='Path to treatment .fastq.gz')
	parser.add_argument('-c', '--control', required=True, type=str, help='Path to control .fastq.gz')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory for processed files')
	parser.add_argument('-n', '--name', required=False, type=str, default='sample', help='Output sample name to prepend')
	parser.add_argument('-ref', '--reference', required=True, type=str, help='Path to reference genome prepared for BWA')
	parser.add_argument('-markdup', '--markdup', required=True, type=str, help='Path to MarkDuplicates.jar')
	parser.add_argument('-p', '--processes', required=False, type=int, default=4, help='Number of processes to use')
	parser.add_argument('-m', '--memory', required=False, type=int, default=8, help='Maximum memory per thread for samtools sort')
	parser.add_argument('-q', '--quality', required=False, type=int, default=10, help='Mapping quality cutoff for samtools')
	parser.add_argument('--broad', required=False, action='store_true', default=False, help='Broad peak option for MACS2 callpeak')
	parser.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip read alignment step')
	parser.add_argument('--skip_peaks', required=False, action='store_true', default=False, help='Skip calling peaks step')
	parser.add_argument('--skip_track', required=False, action='store_true', default=False, help='Skip making signal track for genome browser')
	return parser.parse_args()

#=======================================================#

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
