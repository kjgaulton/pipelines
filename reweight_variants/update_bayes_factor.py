#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

def process_args():
	parser = argparse.ArgumentParser(description='Update PPA for fine-mapping studies using .ridgeparams file')
	parser.add_argument('-r', '--ridgeparams', required=True, type=str, help='Path to ridgeparams output file from fgwas')
	parser.add_argument('-f', '--finemap', required=True, type=str, help='Path to annotated fine-mapping file')
	return parser.parse_args()

def load_ridgeparams(args):
	states = {}
	with open(args.finemap) as f:
		all_states = [x.split('\t')[5] for x in f.read().splitlines()]
		states = {x:0 for x in all_states}
	with open(args.ridgeparams) as f:
		lines = f.read().splitlines()[1:-1]
		for line in lines:
			state, parameter = line.split(' ')
			states[state] = float(parameter)
	return states

def load_finemap(args, states):
	finemap_master = pd.read_table(args.finemap, header=None)
	finemap_master.columns=['chrom', 'pos', 'rsid', 'locus', 'lnbf', 'state']
	finemap_master['parameter'] = [states[x] for x in finemap_master['state']]
	loci = sorted(set(finemap_master['locus']))
	return finemap_master, loci

def update_loci(finemap_master, locus):
	finemap_subset = finemap_master[finemap_master['locus'] == locus]
	finemap_subset['new_lnbf'] = np.log(np.exp(finemap_subset['lnbf']) * (np.exp(finemap_subset['parameter'])/sum(np.exp(finemap_subset['parameter']))))
	finemap_subset['new_ppa'] = np.exp(finemap_subset['new_lnbf'])/sum(np.exp(finemap_subset['new_lnbf']))
	finemap_subset = finemap_subset.sort_values('new_ppa', ascending=False)
	return finemap_subset

def main(args):
	finemap_name = os.path.basename(args.finemap).split('.')[0]
	states = load_ridgeparams(args)
	finemap_master, loci = load_finemap(args, states)
	for locus in loci:
		finemap_subset = update_loci(finemap_master, locus)
		finemap_subset.to_csv('.'.join([locus, finemap_name, 'updated_PPA.tab']), sep='\t', index=False)
	return

args = process_args()
main(args)
