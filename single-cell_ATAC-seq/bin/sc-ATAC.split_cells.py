#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import subprocess
import logging
import traceback


def split_bam(args):
	samtools_view_cmd = ['samtools', 'view', '-h', args.input_bam]
	bam_in = subprocess.Popen(samtools_view_cmd, stdout=subprocess.PIPE)
	header, cells = [], {}
	for line in bam_in.stdout:
		line = line.decode().rstrip('\n')
		if line.startswith(('@HD', '@SQ', '@PG')):
			header.append(line)
		elif line.startswith(('A', 'C', 'G', 'T')):
			barcode = line.split(':')[0]
			cells[barcode] = cells.get(barcode, []) + [line]
		else:
			raise SystemExit('Check if your input bam file is formatted correctly. The bam file needs a header. Each read should start with the barcode.') 
	return header, cells

def filter_barcodes(args, cells):
	barcode_counts = args.prefix + '.barcode_counts.txt'
	barcodes = {b:len(cells[b]) for b in cells.keys()}
	with open(barcode_counts, 'w') as f:
		print('Barcode'.ljust(32), 'Count', sep='\t', file=f)
		for barcode, count in barcodes.items():
			print(barcode, count, sep='\t', file=f)
	filtered = [b for b,c in barcodes.items() if c >= args.min_reads]
	return filtered

def output_cells(args, header, cells, pass_filter):
	reads_before_rmdup = 0
	for bar in pass_filter:
		tmp = tempfile.TemporaryFile('w+b')
		for line in header:
			tmp.write('{}\n'.format(line).encode())
		for line in cells[bar]:
			tmp.write('{}\n'.format(line).encode())
		tmp.seek(0)
		cell_bam = '.'.join([args.prefix, bar, 'filt.bam'])
		samtools_sort_cmd = ['samtools', 'sort', '-m', '{}G'.format(args.memory), '-O', 'bam', '-']
		with open(cell_bam, 'w') as out_bam:
			subprocess.call(samtools_sort_cmd, stdin=tmp, stdout=out_bam)
		tmp.seek(0)
		reads_before_rmdup += sum(1 for _ in tmp)
		tmp.close()
	return reads_before_rmdup / 2 

def remove_duplicates(args, barcodes):
	reads_after_rmdup = 0
	for bar in barcodes:
		input_bam = '.'.join([args.prefix, bar, 'filt.bam'])
		output_bam = '.'.join([args.prefix, bar, 'rmdup.bam'])
		metrics_file = '.'.join([args.prefix, bar, 'picard_rmdup_metrics.txt'])
		rmdup_cmd = [
				'java', '-Xmx{}G'.format(args.memory), 
				'-jar', args.picard_mark_dup,
				'INPUT={}'.format(input_bam),
				'OUTPUT={}'.format(output_bam),
				'REMOVE_DUPLICATES=true',
				'VALIDATION_STRINGENCY=LENIENT',
				'METRICS_FILE={}'.format(metrics_file)
			]
		with open(os.devnull) as f:
			subprocess.call(rmdup_cmd, stderr=f)
		count_reads_cmd = ['samtools', 'view', output_bam]
		p = subprocess.Popen(count_reads_cmd, stdout=subprocess.PIPE)
		reads_after_rmdup += sum(1 for _ in p.stdout)
	return reads_after_rmdup / 2

def main(args):
	logging.info('Starting up.')
	try:
		logging.info('Read in combined bam file [{}]'.format(args.input_bam))
		header, cells = split_bam(args)
		logging.info('Filtering barcodes for [{}] minimum reads.'.format(args.min_reads))
		pass_filter = filter_barcodes(args, cells)
		logging.info('Outputting bams for cells that pass the filter.')
		before_reads = output_cells(args, header, cells, pass_filter)
		logging.info('Removing duplicate reads using picard.')
		after_reads = remove_duplicates(args, pass_filter)
		pass_rate = float(len(pass_filter))/len(cells.keys())
		dup_rate = (before_reads - after_reads) / float(before_reads) * 100
		logging.info('Finishing up.')
		logging.info('Total barcodes: {}'.format(len(cells.keys())))
		logging.info('Barcodes passing filter: {}'.format(len(pass_filter)))
		logging.info('Percentage barcodes passing filter: {0:.2f}%'.format(pass_rate))
		logging.info('Pairs before rmdup: {}'.format(before_reads))
		logging.info('Pairs after rmdup: {}'.format(after_reads))
		logging.info('Est. PCR duplication rate: {0:.2f}%'.format(dup_rate))
	except Exception as e:
		logging.error(e)
		traceback.print_exc(file=sys.stderr)
		sys.exit(1)
	return

def process_args():
	parser = argparse.ArgumentParser(description='Split compiled bam, collect metrics, and remove duplicate read pairs.')
	parser.add_argument('-i', '--input_bam', required=True, type=str, help='Path to combined.bam file')
	parser.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	parser.add_argument('-n', '--name', required=True, type=str, help='Prefix for naming all output files')
	parser.add_argument('-m', '--memory', required=False, type=int, default=4, help='Maximum memory for samtools sort and picard MarkDuplicates.jar')
	parser.add_argument('-min', '--min_reads', required=False, type=int, default=500, help='Minimum number of reads for output cells')
	parser.add_argument('-mark_dup', '--picard_mark_dup', required=False, type=str, default='Picard/MarkDuplicates.jar', help='Path to picard\'s MarkDuplicates.jar tool')
	return parser.parse_args()

args = process_args()
args.prefix = os.path.join(args.output, args.name)
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
