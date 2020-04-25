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
from multiprocessing import Pool

def trim_reads(args):
	pair1_trim = os.path.join(args.output, os.path.basename(args.read1).split('.fastq.gz')[0] + '_val_1.fq.gz')
	pair2_trim = os.path.join(args.output, os.path.basename(args.read2).split('.fastq.gz')[0] + '_val_2.fq.gz')
	if os.path.isfile(pair1_trim) and os.path.isfile(pair2_trim):
		return pair1_trim, pair2_trim
	trim_galore_cmd = ['trim_galore', '--nextera', '--fastqc', '-q', '20', '-o', args.output, '--paired', args.read1, args.read2]
	with open(os.devnull, 'w') as f:
		subprocess.call(trim_galore_cmd, stderr=f, stdout=f)
	return pair1_trim, pair2_trim

def align_reads(args):
	align_log = args.output_prefix + '.align.log'
	aligned_bam = args.output_prefix + '.compiled.bam'
	filtered_bam = args.output_prefix + '.compiled.filt.bam'

	bwa_mem_cmd = ['bwa', 'mem', '-t', str(args.threads), args.reference, args.read1, args.read2]
	filter1_cmd = ['samtools', 'view', '-b', '-h', '-q', str(args.map_quality), '-F', '2048', '-']
	fixmate_cmd = ['samtools', 'fixmate', '-r', '-', '-']
	samtools_sort_cmd = ['samtools', 'sort', '-m', str(args.memory)+'G', '-@', str(args.threads), '-']
	filter2_cmd = ['samtools', 'view', '-h', '-f', '3', '-F', '4', '-F', '256', '-F', '1024', aligned_bam]
	samtools_view_cmd = ['samtools', 'view', '-b', '-']
	
	with open(align_log, 'w') as log, open(aligned_bam, 'w') as bam_out:
		bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE, stderr=log)
		filt1 = subprocess.Popen(filter1_cmd, stdin=bwa_mem.stdout, stdout=subprocess.PIPE)
		fixmate = subprocess.Popen(fixmate_cmd, stdin=filt1.stdout, stdout=subprocess.PIPE)
		subprocess.call(samtools_sort_cmd, stdin=fixmate.stdout, stdout=bam_out, stderr=log)
	
	with open(align_log, 'a') as log, open(filtered_bam, 'w') as bam_out:
		filt2 = subprocess.Popen(filter2_cmd, stderr=log, stdout=subprocess.PIPE)
		view = subprocess.Popen(samtools_view_cmd, stdin=subprocess.PIPE, stdout=bam_out)
		for line in filt2.stdout:
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

def qc_metrics(args):
	md_bam = args.output_prefix + '.filt.md.bam'
	rmdup_bam = args.output_prefix + '.filt.rmdup.bam'
	tagalign_file = args.output_prefix + '.filt.rmdup.tagAlign.gz'
	
	if os.path.isfile(rmdup_bam):
		with gzip.open(tagalign_file, 'wt') as f:
			bamfile = pysam.AlignmentFile(rmdup_bam, 'rb')
			genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
			for read in bamfile:
				if not read.is_proper_pair:
					continue
				barcode = read.query_name.split('_')[0]
				read_chr = read.reference_name
				read_start = max(1, read.reference_end - args.shift - args.extsize - 5 if read.is_reverse else read.reference_start + args.shift + 4)
				read_end = min(genome_size[read_chr], read.reference_end - args.shift - 5 if read.is_reverse else read.reference_start + args.shift + args.extsize + 4)
				read_qual = read.mapping_quality
				if read.is_reverse:
					read_orient = '-'
				else:
					read_orient = '+'
				line_out = '\t'.join([read_chr, str(read_start), str(read_end), '{}_{}'.format(args.name,barcode), str(read_qual), read_orient])
				print(line_out, sep='\t', file=f)
			bamfile.close()
	else:
		raise FileNotFoundError('{} not found!'.format(rmdup_bam))
	
	qc_metrics = {}
	chr_names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
	# mito and duplicated read percentage
	if os.path.isfile(md_bam):
		bamfile = pysam.AlignmentFile(md_bam, 'rb')
		for read in bamfile:
			barcode = '{}_{}'.format(args.name, read.query_name.split('_')[0])
			if barcode not in qc_metrics:
				qc_metrics[barcode] = {}
			if not read.is_duplicate:
				if read.reference_name in chr_names: 
					qc_metrics[barcode]['unique_usable_reads'] = qc_metrics[barcode].get('unique_usable_reads', 0) + 1
				elif read.reference_name == 'chrM':
					qc_metrics[barcode]['unique_mito_reads'] = qc_metrics[barcode].get('unique_mito_reads', 0) + 1
				qc_metrics[barcode]['total_sequenced_reads'] = qc_metrics[barcode].get('total_sequenced_reads', 0) + 1
			else:
				qc_metrics[barcode]['duplicated_reads'] = qc_metrics[barcode].get('duplicated_reads', 0) + 1
				qc_metrics[barcode]['total_sequenced_reads'] = qc_metrics[barcode].get('total_sequenced_reads', 0) + 1
	else:
		raise FileNotFoundError('{} not found!'.format(md_bam))
	qc_metrics = pd.DataFrame.from_dict(qc_metrics, orient='index').fillna(0).astype(int)
	# reads in peaks
	macs2_cmd = ['macs2', 'callpeak', '-t', tagalign_file, '--outdir', args.output, '-n', args.name, '-p', '.05', '--nomodel', '--keep-dup', 'all', '--shift', '0', '--extsize', '200', '-g', 'hs']
	with open(os.devnull, 'w') as f:
		subprocess.call(macs2_cmd, stderr=f)
	try:
		os.remove(args.output_prefix + '_peaks.xls')
		os.remove(args.output_prefix + '_summits.bed')
	except:
		pass
	blacklist_cmd = subprocess.Popen(['bedtools', 'intersect', '-a' , args.output_prefix + '_peaks.narrowPeak', '-b', args.blacklist_file, '-v'], stdout=subprocess.PIPE)
	intersect_cmd = subprocess.Popen(['bedtools', 'intersect', '-a', tagalign_file, '-b', '-' ], stdin=blacklist_cmd.stdout, stdout=subprocess.PIPE)
	peak_counts = {bc:0 for bc in qc_metrics.index}
	for line in intersect_cmd.stdout:
		line = line.decode()
		fields = line.rstrip().split('\t')
		peak_counts[fields[3]] += 1
	qc_metrics['reads_in_peaks'] = qc_metrics.index.map(peak_counts).fillna(0).astype(int)
	# reads in promoter and promoter usage
	tss_counts = {bc:0 for bc in qc_metrics.index}
	tss_cmd1 = subprocess.Popen(['bedtools', 'intersect', '-a', tagalign_file, '-b', args.promoter_file, '-wa'], stdout=subprocess.PIPE)
	uniq = subprocess.Popen(['uniq'], stdin=tss_cmd1.stdout, stdout=subprocess.PIPE)
	for line in uniq.stdout:
		line = line.decode()
		fields = line.rstrip().split('\t')
		tss_counts[fields[3]] += 1

	tss_used = {bc:[] for bc in qc_metrics.index}
	tss_cmd2 = subprocess.Popen(['bedtools', 'intersect', '-a', tagalign_file, '-b', args.promoter_file, '-wa', '-wb'], stdout=subprocess.PIPE)
	for line in tss_cmd2.stdout:
		line = line.decode()
		fields = line.rstrip().split('\t')
		tss_used[fields[3]].append(fields[9])
	qc_metrics['reads_in_promoters'] = qc_metrics.index.map(tss_counts).fillna(0).astype(int)
	qc_metrics['tss_used'] = [len(set(tss_used[bc])) for bc in qc_metrics.index]
	total_prom = len(sorted(set(pd.read_table(args.promoter_file, sep='\t', header=None, names=['chr','start','end','promoter'])['promoter'])))
	qc_metrics['frac_reads_in_peaks'] = qc_metrics['reads_in_peaks'].div(qc_metrics['unique_usable_reads']).replace(np.inf, 0).fillna(0)
	qc_metrics['frac_reads_in_promoters'] = qc_metrics['reads_in_promoters'].div(qc_metrics['unique_usable_reads']).replace(np.inf, 0).fillna(0)
	qc_metrics['frac_promoters_used'] = qc_metrics['tss_used']/total_prom
	qc_metrics['frac_mito_reads'] = qc_metrics['unique_mito_reads'].div(qc_metrics['unique_usable_reads'] + qc_metrics['unique_mito_reads']).replace(np.inf, 0).fillna(0)
	qc_metrics['frac_duplicated_reads'] = qc_metrics['duplicated_reads'].div(qc_metrics['total_sequenced_reads']).fillna(0)
	qc_metrics.to_csv(os.path.join(args.output_prefix + '.qc_metrics.txt'), sep='\t')
	return	


def generate_matrix(args):
	tagalign_file = args.output_prefix + '.filt.rmdup.tagAlign.gz'
	qc_metrics = pd.read_table(args.output_prefix + '.qc_metrics.txt', sep='\t', header=0, index_col=0)
	pass_barcodes = qc_metrics.loc[qc_metrics['unique_usable_reads'] >= args.minimum_reads].index 

	lf_mtx_file = args.output_prefix + '.long_fmt_mtx.txt.gz'
	barcodes_file = args.output_prefix + '.barcodes'
	regions_file = args.output_prefix + '.regions'
	mtx_file = args.output_prefix + '.mtx'
	
	windows_file = make_windows(args)
	window_intersect = intersect_regions(tagalign_file, windows_file)
	cut = subprocess.Popen(['cut', '-f', '4,10'], stdin=window_intersect.stdout, stdout=subprocess.PIPE)
	sort = subprocess.Popen(['sort', '-S', '{}G'.format(args.memory * args.threads)], stdin=cut.stdout, stdout=subprocess.PIPE)
	uniq = subprocess.Popen(['uniq', '-c'], stdin=sort.stdout, stdout=subprocess.PIPE)
	awk = subprocess.Popen(['awk', '''BEGIN{{OFS="\\t"}} {{print $2,$3,$1}}'''], stdin=uniq.stdout, stdout=subprocess.PIPE)
	with gzip.open(lf_mtx_file, 'wt') as f:
		subprocess.call(['gzip', '-c'], stdin=awk.stdout, stdout=f)
	
	lf_mtx = pd.read_table(lf_mtx_file, sep='\t', header=None, names=['barcode','region','count'])
	lf_mtx = lf_mtx.loc[lf_mtx['barcode'].isin(pass_barcodes)]
	lf_mtx.to_csv(lf_mtx_file, sep='\t', header=False, index=False, compression='gzip')

	tmp_R = args.output_prefix + '.tmp.R'
	with open(tmp_R, 'w') as tR:
		print('''library(Matrix)''', file=tR)
		print('''data <- read.table('{}', sep='\\t', header=FALSE)'''.format(lf_mtx_file), file=tR)
		print('''sparse.data <- with(data, sparseMatrix(i=as.numeric(V1), j=as.numeric(V2), x=V3, dimnames=list(levels(V1), levels(V2))))''', file=tR)
		print('''t <- writeMM(sparse.data, '{}')'''.format(mtx_file), file=tR)
		print('''write.table(data.frame(rownames(sparse.data)), file='{}', col.names=FALSE, row.names=FALSE, quote=FALSE)'''.format(barcodes_file), file=tR)
		print('''write.table(data.frame(colnames(sparse.data)), file='{}', col.names=FALSE, row.names=FALSE, quote=FALSE)'''.format(regions_file), file=tR)
	subprocess.call(['Rscript', tmp_R])
	subprocess.call(['gzip', mtx_file])
	os.remove(windows_file)
	os.remove(tmp_R)
	return


def intersect_regions(tagalign, regions):
	awk_cmd = ['awk', '''BEGIN{{FS=OFS="\\t"}} {{peakid=$1":"$2"-"$3; gsub("chr","",peakid); print $1, $2, $3, peakid}}''', regions]
	intersect_cmd = ['bedtools', 'intersect', '-a', tagalign, '-b', '-', '-wa', '-wb']
	awk = subprocess.Popen(awk_cmd, stdout=subprocess.PIPE)
	intersect = subprocess.Popen(intersect_cmd, stdin=awk.stdout, stdout=subprocess.PIPE)
	return intersect

def make_windows(args):
	makewindows_cmd = ['bedtools', 'makewindows', '-g', args.chrom_sizes, '-w', str(args.window_size * 1000)]
	filter_cmd = ['grep', '-v', '_']
	blacklist_cmd = ['bedtools', 'intersect', '-a', '-', '-b', args.blacklist_file, '-v']
	windows_file = '{}.{}kb_windows.bed'.format(args.output_prefix, args.window_size)
	with open(windows_file, 'w') as f:
		makewindows = subprocess.Popen(makewindows_cmd, stdout=subprocess.PIPE)
		filt = subprocess.Popen(filter_cmd, stdin=makewindows.stdout, stdout=subprocess.PIPE)
		subprocess.call(blacklist_cmd, stdin=filt.stdout, stdout=f)
	return windows_file

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
	if not args.skip_qc:
		qc_metrics(args)
	if not args.skip_matrix:
		logging.info('Generating tagalign and chromatin accessibility matrix.')
		generate_matrix(args)
	logging.info('Finish.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Process combinatorial barcoding snATAC-seq data.')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-r1', '--read1', required=True, type=str, help='Paired-end reads file 1')
	io_group.add_argument('-r2', '--read2', required=True, type=str, help='Paired-end reads file 2')
	io_group.add_argument('-o', '--output', required=False, type=str, default=os.getcwd(), help='Output directory to store processed files')
	io_group.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	
	align_group = parser.add_argument_group('Alignment arguments')
	align_group.add_argument('-t', '--threads', required=False, type=int, default=8, help='Number of threads to use for alignment [8]')
	align_group.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum amount of memory (G) per thread for samtools sort [4]')
	align_group.add_argument('-q', '--map_quality', required=False, type=int, default=30, help='Mapping quality score filter for samtools [30]')
	align_group.add_argument('-ref', '--reference', required=False, type=str, default='references/male.hg19.fa', help='Path to the BWA indexed reference genome')
	
	dup_group = parser.add_argument_group('Remove duplicates arguments')
	dup_group.add_argument('--picard', required=False, type=str, default='/home/joshchiou/bin/picard.jar', help='Path to picard.jar')
	
	matrix_group = parser.add_argument_group('Matrix generation arguments')
	matrix_group.add_argument('--shift', required=False, type=int, default=-100, help='Read shift size')
	matrix_group.add_argument('--extsize', required=False, type=int, default=200, help='Read extension size')
	matrix_group.add_argument('--minimum-reads', required=False, type=int, default=1000, help='Minimum number of reads for barcode inclusion')
	matrix_group.add_argument('--window-size', required=False, type=int, default=5, help='Size (kb) to use for defining windows of accessibility')
	matrix_group.add_argument('--chrom-sizes', required=False, type=str, default='references/hg19.chrom.sizes', help='Chromosome sizes file from UCSC')
	matrix_group.add_argument('--blacklist-file', required=False, type=str, default='references/hg19-blacklist.v1.bed.gz', help='BED file of blacklisted regions')
	matrix_group.add_argument('--promoter-file', required=False, type=str, default='references/gencode.v19.1kb_promoters.pc_transcripts.bed.gz', help='BED file of autosomal promoter regions')

	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-trim', required=False, action='store_true', default=False, help='Skip adapter trimming step')
	skip_group.add_argument('--skip-align', required=False, action='store_true', default=False, help='Skip read alignment step')
	skip_group.add_argument('--skip-rmdup', required=False, action='store_true', default=False, help='Skip duplicate removal step')
	skip_group.add_argument('--skip-qc', required=False, action='store_true', default=False, help='Skip QC metrics calculation step')
	skip_group.add_argument('--skip-matrix', required=False, action='store_true', default=False, help='Skip matrix generation step')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
	args = process_args()
	main(args)
