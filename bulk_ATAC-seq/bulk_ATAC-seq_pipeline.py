#!/usr/bin/env python3

import argparse
import subprocess
import os 
import sys
import reprlib
import logging

#=======================================================#

def detect_file_format(args):
	if args.paired1.endswith('.gz') and args.paired2.endswith('.gz'):
		return args.paired1, args.paired2
	if args.paired1.endswith('.bz2') and args.paired2.endswith('.bz2'):
		logging.info('BZ2 file format detected -- converting to GZ file format.')
		p1_bn, p2_bn = args.paired1.split('.bz2')[0], args.paired2.split('.bz2')[0]
		subprocess.call(['bunzip2', args.paired1, args.paired2])
		subprocess.call(['gzip', p1_bn, p2_bn])
		args.paired1, args.paired2 = p1_bn + '.gz', p2_bn + '.gz'
		return args.paired1, args.paired2
	if args.paired1.endswith('.fastq') and args.paired2.endswith('.fastq'):
		logging.info('Unzipped FASTQ format detected -- converting to GZ file format.')
		subprocess.call(['gzip', args.paired1, args.paired2])
		args.paired1, args.paired2 = args.paired1 + '.gz', args.paired2 + '.gz'
		return args.paired1, args.paired2
	if args.paired1.endswith('.fq') and args.paired2.endswith('.fq'):
		logging.info('Unzipped FQ format detected -- converting to GZ file format.')
		p1_bn, p2_bn = args.paired1.split('.fq')[0] + '.fastq', args.paired2.split('.fq')[0] + '.fastq'
		os.rename(args.paired1, p1_bn)
		os.rename(args.paired2, p2_bn)
		subprocess.call(['gzip', p1_bn, p2_bn])
		args.paired1, args.paired2 = p1_bn + '.gz', p2_bn + '.gz'
		return args.paired1, args.paired2
	else:
		logging.error('Unknown file format or paired end reads have different file formats.')
		raise RuntimeError('Paired end reads must be in a valid file format!')
	return

#=======================================================#

def trim_galore(args):
	trim_galore_cmd = ['trim_galore', '--fastqc', '-q', '10', '-o', args.output, '--paired', args.paired1, args.paired2]
	with open(os.devnull, 'w') as f:
		subprocess.call(trim_galore_cmd, stderr=f)
	trim_output_1 = os.path.join(args.output, os.path.basename(args.paired1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	trim_output_2 = os.path.join(args.output, os.path.basename(args.paired2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	return trim_output_1, trim_output_2

#=======================================================#

def process_reads(args):
	output_prefix = os.path.join(args.output, args.name)
	align_log = output_prefix + '.align.log'
	aligned_bam = output_prefix + '.sort.filt.bam'
	rmdup_bam = output_prefix + '.sort.filt.rmdup.bam'

	bwa_mem_cmd = ['bwa', 'mem', '-M', '-t', str(args.threads), args.reference, args.paired1, args.paired2]
	quality_filt_cmd = ['samtools', 'view', '-h', '-f', '0x2', '-F', '0x100', '-q', str(args.quality), '-@', str(args.threads), '-']
	mito_filt_cmd = ['grep', '-v', 'chrM']
	samtools_sort_cmd = ['samtools', 'sort', '-m', '{}G'.format(args.memory), '-@', str(args.threads), '-']
	
	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		qual_filt = subprocess.Popen(quality_filt_cmd, stdin=bwa_mem.stdout, stdout=subprocess.PIPE, stderr=log)
		mito_filt = subprocess.Popen(mito_filt_cmd, stdin=qual_filt.stdout, stdout=subprocess.PIPE, stderr=log)
		subprocess.call(samtools_sort_cmd, stdin=mito_filt.stdout, stdout=bam_out, stderr=log)
	
	
	metrics_file = output_prefix + '.picard_rmdup_metrics.txt'
	rmdup_cmd = [
			'java', '-Xmx{}G'.format(args.memory), 
			'-jar', args.picard_mark_dup,
			'INPUT={}'.format(aligned_bam),
			'OUTPUT={}'.format(rmdup_bam),
			'REMOVE_DUPLICATES=true',
			'VALIDATION_STRINGENCY=LENIENT',
			'METRICS_FILE={}'.format(metrics_file)
		]
	index_cmd = ['samtools', 'index', rmdup_bam]
	if os.path.exists(aligned_bam) and os.path.getsize(aligned_bam) != 0:
		with open(os.devnull, 'w') as f:
			subprocess.call(rmdup_cmd, stderr=f)
	if os.path.exists(rmdup_bam) and os.path.getsize(rmdup_bam) != 0:
		subprocess.call(index_cmd)
	
	return rmdup_bam

#=======================================================#

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

def make_bedgraph(args):
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
	
	label = 'track type=bedGraph name=\"{0}\" description=\"{0}\" visibility=2 color=0,0,0 altColor=0,0,0 autoScale=off maxHeightPixels=64:64:32'.format(args.name)
	sorted_bdg = os.path.join(args.output, args.name  + '_ppois.sorted.bdg')
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
	if not os.path.isdir(args.output):
		try:
			os.makedirs(args.output)
		except OSError:
			pass
	
	args.paired1, args.paired2 = detect_file_format(args)
	
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
	if not args.skip_bdg:
		try:
			logging.info('Generating signal track with MACS2')
			make_bedgraph(args)
		except:
			print(repr(e), file=sys.stderr)
	logging.info('Cleaning up temporary files.')
	if not args.skip_cleanup:
		cleanup(args)
	logging.info('Finishing up.')
	return

#=======================================================#

def process_args():
	parser = argparse.ArgumentParser(description='Pipeline for ATAC-seq to trim reads, align them to a reference genome, call peaks, and generate a genome browser signal track.')
	
	io_group = parser.add_argument_group('I/O arguments')	
	io_group.add_argument('-p1', '--paired1', required=False, type=str, help='Path to paired reads (1)')
	io_group.add_argument('-p2', '--paired2', required=False, type=str, help='Path to paired reads (2)')
	io_group.add_argument('-o', '--output', required=True, type=str, help='Output directory for processed files')
	io_group.add_argument('-n', '--name', required=False, type=str, default='sample', help='Output sample name to prepend')
	
	align_group = parser.add_argument_group('Alignment arguments')
	align_group.add_argument('-t', '--threads', required=False, type=int, default=4, help='Number of threads to use [4]')
	align_group.add_argument('-m', '--memory', required=False, type=int, default=8, help='Maximum memory per thread for samtools sort [8]')
	align_group.add_argument('-q', '--quality', required=False, type=int, default=30, help='Mapping quality cutoff for samtools [30]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='/home/joshchiou/references/ucsc.hg19.fasta', help='Path to reference genome [/home/joshchiou/references/ucsc.hg19.fasta]')
	align_group.add_argument('--picard_mark_dup', required=False, type=str, default='/home/joshchiou/bin/MarkDuplicates.jar', help='Path to picard MarkDuplicates.jar [/home/joshchiou/bin/MarkDuplicates.jar]')

	skip_group = parser.add_argument_group('Skip processing steps')
	skip_group.add_argument('--skip_trim', required=False, action='store_true', default=False, help='Skip adapter trimming step')
	skip_group.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip read alignment step')
	skip_group.add_argument('--skip_peaks', required=False, action='store_true', default=False, help='Skip calling peaks step')
	skip_group.add_argument('--skip_bdg', required=False, action='store_true', default=False, help='Skip making genome browser track')
	skip_group.add_argument('--skip_cleanup', required=False, action='store_true', default=False, help='Skip cleanup operations')
	
	return parser.parse_args()

#=======================================================#
if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
