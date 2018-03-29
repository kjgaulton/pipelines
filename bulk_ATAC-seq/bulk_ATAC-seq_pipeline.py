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
		subprocess.call(trim_galore_cmd, stderr=f, stdout=f)
	trim_output_1 = os.path.join(args.output, os.path.basename(args.paired1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	trim_output_2 = os.path.join(args.output, os.path.basename(args.paired2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	return trim_output_1, trim_output_2

#=======================================================#

def process_reads(args):
	output_prefix = os.path.join(args.output, args.name)
	align_log = output_prefix + '.align.log'
	sort_bam = output_prefix + '.sort.bam'
	md_bam = output_prefix + '.sort.md.bam'
	rmdup_bam = output_prefix + '.sort.filt.rmdup.bam'

	bwa_mem_cmd = ['bwa', 'mem', '-M', '-t', str(args.threads), args.reference, args.paired1, args.paired2]
	samtools_sort_cmd = ['samtools', 'sort', '-m', '{}G'.format(args.memory), '-@', str(args.threads), '-']

	with open(align_log, 'w') as log, open(sort_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		subprocess.call(samtools_sort_cmd, stdin=bwa_mem.stdout, stdout=bam_out, stderr=log)
	
	metrics_file = output_prefix + '.picard_rmdup_metrics.txt'
	
	md_cmd = [
			'java', '-Xmx{}G'.format(args.memory), 
			'-jar', args.picard_mark_dup,
			'INPUT={}'.format(sort_bam),
			'OUTPUT={}'.format(md_bam),
			'REMOVE_DUPLICATES=false',
			'ASSUME_SORTED=true',
			'VALIDATION_STRINGENCY=LENIENT',
			'METRICS_FILE={}'.format(metrics_file)]
	rmdup_cmd = [
			'samtools', 'view',
			'-b', '-h', '-f', '3',
			'-F', '4', '-F', '256',
			'-F', '1024', '-F', '2048',
			'-q', str(args.quality),
			md_bam]

	autosomal_chr = ['chr' + str(c) for c in range(1,23)]
	rmdup_cmd.extend(autosomal_chr)

	if os.path.exists(sort_bam) and os.path.getsize(sort_bam) != 0:
		with open(os.devnull, 'w') as f:
			subprocess.call(md_cmd, stderr=f)
	if os.path.exists(md_bam) and os.path.getsize(md_bam) != 0:
		subprocess.call(['samtools', 'index', md_bam])
		with open(rmdup_bam, 'w') as f:
			subprocess.call(rmdup_cmd, stdout=f)
		subprocess.call(['samtools', 'index', rmdup_bam])

	return md_bam, rmdup_bam

#=======================================================#

def call_peaks(args, input_bam):
	macs2_log = os.path.join(args.output, '.'.join([args.name, 'macs2_callpeaks.log']))
	macs2_cmd = ['macs2', 'callpeak', 
			'-t', input_bam, 
			'--outdir', args.output, 
			'-n', args.name,
			'-g', args.macs2_genome,
			'--nomodel', '-B', 
			'--shift', '-100', 
			'--extsize', '200', 
			'--keep-dup', 'all']
	with open(macs2_log, 'w') as f:
		subprocess.call(macs2_cmd, stderr=f)
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
	
	bdg_label = 'track type=bedGraph name=\"{0}\" description=\"{0}, ATAC signal\" visibility=2 color={1} altColor=0,0,0 autoScale=off maxHeightPixels=64:64:32'.format(args.name, args.bdg_color)
	sorted_bdg = os.path.join(args.output, args.name  + '_ppois.sorted.bdg')
	sort_cmd = ['sort', '-k', '1,1', '-k', '2,2n', bdgcmp_out]
	with open(sorted_bdg, 'w') as f:
		print(bdg_label, file=f)
	with open(sorted_bdg, 'a') as f:	
		subprocess.call(sort_cmd, stdout=f)
	subprocess.call(['gzip', sorted_bdg])
	os.remove(bdgcmp_out)
	
	narrow_peaks = os.path.join(args.output, args.name + '_peaks.narrowPeak')
	peak_label = 'track type=narrowPeak name=\"{0} peaks\" description=\"{0}, ATAC peaks\" color={1}'.format(args.name, args.bdg_color)
	peaks = open(narrow_peaks).read().splitlines()
	with open(narrow_peaks, 'w') as f:
		print(peak_label, file=f)
		for p in peaks:
			print(p, file=f)

	return

#=======================================================#

def get_qc_metrics(args, md_bam):
	ATAC_peaks = os.path.join(args.output, args.name + '_peaks.narrowPeak')
	qc_json = os.path.join(args.output, args.name + '.ataqv.json.gz')
	qc_log = os.path.join(args.output, args.name + '.ataqv.log')
	qc_cmd = [
			'ataqv', '--verbose',
			'--metrics-file', qc_json,
			'--peak-file', ATAC_peaks,
			'--tss-file', args.tss,
			'--tss-extension', '1000',
			'--excluded-region-file', args.blacklist,
			'--name', args.name,
			'--description', 'Gaulton lab ATAC sample {}'.format(args.name),
			'--mitochondrial-reference-name', 'chrM',
			'--threads', str(args.threads),
			'human', md_bam]
	with open(qc_log, 'w') as f, open(os.devnull, 'w') as n:
		subprocess.call(qc_cmd, stdout=f, stderr=n)
	return

#=======================================================#

def cleanup(args):
	read_1_name = os.path.basename(args.paired1).split('.fastq.gz')[0]
	read_2_name = os.path.basename(args.paired2).split('.fastq.gz')[0]
	#trim_tmp_1 = os.path.join(args.output, read_1_name + '_trimmed.fq.gz')
	#trim_tmp_2 = os.path.join(args.output, read_2_name + '_trimmed.fq.gz')
	fastqc_zip_1 = os.path.join(args.output, read_1_name + '_val_1_fastqc.zip')
	fastqc_zip_2 = os.path.join(args.output, read_2_name + '_val_2_fastqc.zip')
	clean = [fastqc_zip_1, fastqc_zip_2]
	try:
		for f in clean:
			os.remove(f)
	except:
		pass
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
			logging.error('Failed during adaptor trimming step with trim_galore.')
			print('Check options -p1 and -p2: ' + repr(e), file=sys.stderr)
			sys.exit(1)
	if not args.skip_align:
		try:
			logging.info('Aligning reads with bwa and filtering reads with samtools.')
			md_bam, rmdup_bam = process_reads(args)
		except Exception as e:
			logging.error('Failed during alignment step with bwa mem.')
			print(repr(e), file=sys.stderr)
			sys.exit(1)
	if not args.skip_peaks:
		try:
			logging.info('Calling peaks with MACS2.')
			call_peaks(args, rmdup_bam)
		except Exception as e:
			logging.error('Failed during MACS2 peak calling step.')
			print(repr(e), file=sys.stderr)
			sys.exit(1)
	if not args.skip_bdg:
		try:
			logging.info('Generating signal track with MACS2.')
			make_bedgraph(args)
		except Exception as e:
			logging.error('Failed during bedgraph generation step.')
			print(repr(e), file=sys.stderr)
			sys.exit(1)
	logging.info('Cleaning up temporary files.')
	if not args.skip_qc:
		try:
			logging.info('Starting up QC using ataqv.')
			get_qc_metrics(args, md_bam)
		except Exception as e:
			logging.error('Failed during QC step.')
			print(repr(e), file=sys.stderr)
			sys.exit(1)
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
	align_group.add_argument('-m', '--memory', required=False, type=int, default=8, help='Maximum memory (in Gb) per thread for samtools sort [8]')
	align_group.add_argument('-q', '--quality', required=False, type=int, default=30, help='Mapping quality cutoff for samtools [30]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='/home/joshchiou/references/male.hg19.fa', help='Path to reference genome [/home/joshchiou/references/male.hg19.fa]')
	align_group.add_argument('--picard_mark_dup', required=False, type=str, default='/home/joshchiou/bin/MarkDuplicates.jar', help='Path to picard MarkDuplicates.jar [/home/joshchiou/bin/MarkDuplicates.jar]')
	
	qc_group = parser.add_argument_group('QC arguments')
	qc_group.add_argument('--tss', required=False, type=str, default='/home/joshchiou/references/hg19_gencode_tss_unique.bed', help='Path to TSS definitions for calculating ATAC signal enrichment around TSS [/home/joshchiou/references/hg19_gencode_tss_unique.bed]')
	qc_group.add_argument('--blacklist', required=False, type=str, default='/home/joshchiou/references/ENCODE.hg19.blacklist.bed', help='Path to blacklist BED file to ignore ENCODE high signal regions [/home/joshchiou/references/ENCODE.hg19.blacklist.bed]')
	
	bedgraph_group = parser.add_argument_group('Signal track arguments')
	bedgraph_group.add_argument('--macs2_genome', required=False, type=str, default='hs', help='MACS2 genome (e.g. hs or mm) for peak calling')
	bedgraph_group.add_argument('--bdg_color', required=False, type=str, default='0,0,0', help='Color for genome browser signal track in R,G,B [0,0,0]')

	skip_group = parser.add_argument_group('Skip processing steps')
	skip_group.add_argument('--skip_trim', required=False, action='store_true', default=False, help='Skip adapter trimming step [OFF]')
	skip_group.add_argument('--skip_align', required=False, action='store_true', default=False, help='Skip read alignment step [OFF]')
	skip_group.add_argument('--skip_peaks', required=False, action='store_true', default=False, help='Skip calling peaks step [OFF]')
	skip_group.add_argument('--skip_bdg', required=False, action='store_true', default=False, help='Skip making genome browser track [OFF]')
	skip_group.add_argument('--skip_qc', required=False, action='store_true', default=False, help='Skip ATAC qc step using ataqv [OFF]')
	skip_group.add_argument('--skip_cleanup', required=False, action='store_true', default=False, help='Skip cleanup operations [OFF]')
	
	return parser.parse_args()

#=======================================================#
if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
