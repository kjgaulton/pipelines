#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
pd.options.display.float_format = '${:,.2f}'.format

def process_args():
	parser = argparse.ArgumentParser(description='Insert a useful description here')
	parser.add_argument('-d','--dirs', nargs='+', required=True, type=str, help='List of directories containing footprints')
	parser.add_argument('-l', '--motifs', required=False, type=str, default='/home/joshchiou/references/combined_motifs.list', help='Path to list of motifs')
	return parser.parse_args()

def proc_header(args):
	basenames = [os.path.basename(d.rstrip('/')) for d in args.dirs]
	return basenames

def load_motifs(args):
	with open(args.motifs) as f:
		motif_list = f.read().splitlines()
	return motif_list

def count_footprints(fp_file):
	if not os.path.exists(fp_file):
		return 0
	with open(fp_file) as f:
		count = sum(1 for line in f)
	return count

def process_dataframe(table):
	df = pd.DataFrame(table[1:], columns=table[0]).set_index('Motif')
	df = df[(df.T != 0).any()]
	totals = [np.sum(df[column]) for column in df]
	norm_factor = totals/min(totals)
	for i,column in enumerate(df):
		df[column] = df[column] / norm_factor[i]
	return df

def main(args):
	table = [['Motif'] + proc_header(args)]
	motif_list = load_motifs(args)
	for motif in motif_list:
		counts = [count_footprints(os.path.join(d, '.'.join([motif,os.path.basename(d.rstrip('/')),'footprints.bed']))) for d in args.dirs]
		table.append([motif] + counts)
	data = process_dataframe(table)
	data.round(2).to_csv('normalized_counts.tab', sep='\t')
	return

args = process_args()
main(args)
