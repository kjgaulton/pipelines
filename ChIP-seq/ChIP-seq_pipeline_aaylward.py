#!/usr/bin/env python3

import argparse
import subprocess
import os 
import sys
import reprlib
import logging
import socket

hostname = socket.gethostname()
if hostname == 'gatsby.ucsd.edu':
    sys.path.append('/home/data/kglab-python3-modules')
elif hostname == 'holden':
    sys.path.append('/lab/kglab-python3-modules')

import seqalign

#=======================================================#

def process_reads(args, reads, name):
	output_prefix = os.path.join(args.output, name)
	align_log = output_prefix + '.align.log'
	aligned_bam = output_prefix + '.sort.filt.bam'
	rmdup_bam = output_prefix + '.sort.filt.rmdup.bam'
	
	# align reads with bwa aln / bwa samse
	# then filter out reads that are unmapped, < quality score
	# then sort reads
	# then remove chrM reads
	with open(align_log, 'w') as log:
		with seqalign.SequenceAlignment(
			input_file=reads,
			phred_quality_score=args.quality,
			processes=args.processes,
			log=log,
			aligner=seqalign.BwaAligner(trim=15)
		) as sa:
			sa.cleans_up_bam=False
			sa.apply_quality_filter()
			sa.samtools_sort(memory_limit=args.memory * args.processes)
			sa.samtools_index()
			sa.remove_mitochondrial_reads()
			sa.write(aligned_bam)
	
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
			'-q', str(args.qvalue),
			'--extsize', '200', 
			'-B', 
			'--keep-dup', 'all' ]
	if args.broad:
		macs2_cmd.extend(['--broad', '--broad-cutoff', str(args.broad_cutoff)])
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
	label = 'track type=bedGraph name=\"{0}\" description=\"{0}\" visibility=2 color={1} altColor=0,0,0 autoScale=off maxHeightPixels=64:64:32'.format(args.name, args.color)
	sorted_bdg = os.path.join(args.output, args.name  + '_ppois.sorted.bdg')
	sort_cmd = ['sort', '-k', '1,1', '-k', '2,2n', bdgcmp_out]
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
	if not os.path.isdir(args.output):
		try:
			os.makedirs(args.output)
		except OSError:
			pass
	if not args.skip_align:
		try:
			logging.info('Aligning reads with bwa aln/samse and filtering reads with samtools.')
			logging.info('Processing treatment bam [{}].'.format(args.treatment))
			args.treatment = process_reads(args, args.treatment, args.name)
			logging.info('Processing control bam [{}].'.format(args.control))
			args.control = process_reads(args, args.control, args.name + '_control')
		except Exception as e:
			print('Check options -t and -c: ' + repr(e), file=sys.stderr)
	if not args.skip_peaks:
		try:
			logging.info('Calling peaks with MACS2.')
			call_peaks(args, args.treatment, args.control)
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
	
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-t', '--treatment', required=True, type=str, help='Path to treatment file [.fastq.gz OR .bam if --skip_align is ON]')
	io_group.add_argument('-c', '--control', required=True, type=str, help='Path to control file [.fastq.gz OR .bam if --skip_align is ON]')
	io_group.add_argument('-o', '--output', required=True, type=str, help='Output directory for processed files')
	io_group.add_argument('-n', '--name', required=False, type=str, default='sample', help='Output sample name to prepend')
	
	align_group = parser.add_argument_group('Alignment and rmdup arguments')
	align_group.add_argument('-p', '--processes', required=False, type=int, default=4, help='Number of processes to use [4]')
	align_group.add_argument('-m', '--memory', required=False, type=int, default=8, help='Maximum memory per thread [8]')
	align_group.add_argument('-q', '--quality', required=False, type=int, default=10, help='Mapping quality cutoff for samtools [10]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='/home/joshchiou/references/ucsc.hg19.fasta',  help='Path to reference genome prepared for BWA [/home/joshchiou/references/ucsc.hg19.fasta]')
	align_group.add_argument('-markdup', '--markdup', required=False, type=str, default='/home/joshchiou/bin/MarkDuplicates.jar', help='Path to MarkDuplicates.jar [/home/joshchiou/bin/MarkDuplicates.jar]')
	
	macs2_group = parser.add_argument_group('MACS2 parameters')
	macs2_group.add_argument('--qvalue', required=False, type=float, default=0.01, help='MACS2 callpeak qvalue cutoff [0.01]')
	macs2_group.add_argument('--broad', required=False, action='store_true', default=False, help='Broad peak option for MACS2 callpeak [OFF]')
	macs2_group.add_argument('--broad_cutoff', required=False, type=float, default=0.05, help='MACS2 callpeak qvalue cutoff for broad regions [0.05]')
	macs2_group.add_argument('--color', required=False, type=str, default='0,0,0', help='Color in R,G,B format to display for genome browser track [0,0,0]')

	skip_group = parser.add_argument_group('Skip processing')
	skip_group.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip read alignment step [OFF]')
	skip_group.add_argument('--skip_peaks', required=False, action='store_true', default=False, help='Skip calling peaks step [OFF]')
	skip_group.add_argument('--skip_track', required=False, action='store_true', default=False, help='Skip making signal track for genome browser [OFF]')
	return parser.parse_args()

#=======================================================#
if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
