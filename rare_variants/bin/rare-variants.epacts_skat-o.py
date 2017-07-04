#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import logging
from multiprocessing import Pool


def run_epacts(args, name):
    vcf_file = os.path.join(args.proc_vcf_dir, name + '.vcf.gz')
    group_file = os.path.join(args.group_dir, name + '.group')
    out_prefix = os.path.join(args.out_dir, name)
    epacts_cmd = [
            'epacts', 'group',
            '--vcf', vcf_file,
            '--groupf', group_file,
            '--out', out_prefix,
            '--ped', args.ped_file,
            '--test', 'skat', '--skat-o',
            '--max-MAF', str(args.maximum_maf),
            '--min-rsq', str(args.minimum_rsquared),
            '--pheno', 'T2D',
            '--cov', 'Age',
            '--cov', 'Sex',
            '--cov', 'PC_GoT2D_1',
            '--cov', 'PC_GoT2D_2',
            '--cov', 'PC_GoT2D_3',
            '--cov', 'NUMSING',
            '--cov', 'NUMRARE']
    try:
        with open(os.devnull, 'w') as f:
            subprocess.call(epacts_cmd, stderr=f)
            subprocess.call(['make', '-f', out_prefix + '.Makefile', '-j', '1'], stderr=f)
            for run_file in [out_prefix+'.phe', out_prefix+'.ind', out_prefix+'.Makefile', out_prefix+'.cov', out_prefix+'.epacts.OK']:
                os.remove(run_file)
    except Exception as e:
        logging.error(e)
        pass
    return

def main(args):
    logging.info('Starting up.')
    with open('fp.list') as f:
        fp_names = f.read().splitlines()
    args_list = list(zip([args for i in fp_names], fp_names))
    pool = Pool(processes=(args.threads))
    pool.starmap(process_vcf, args_list)
    pool.close()
    pool.join()
    logging.info('Finishing up.')
    return

def process_args():
    parser = argparse.ArgumentParser(description='Run epacts for each group.')
    parser.add_argument('--proc_vcf_dir', required=False, type=str, default='proc_vcfs/', help='Path to processed variants directory')
    parser.add_argument('--group_dir', required=False, type=str, default='groups/', help='Path to groups directory containing group files')
    parser.add_argument('--out_dir', required=False, type=str, default='skat-o/', help='Path to out directory for epacts group test')
    parser.add_argument('--ped_file', required=False, type=str, default='/home/joshchiou/T2D/GoT2D.PC.ped', help='Path to PED file containing covariates')
    parser.add_argument('-min_rsq', '--minimum_rsquared', required=False, type=float, default=0.6, help='Minimum R^2 value to filter variants')
    parser.add_argument('-max_maf', '--maximum_maf', required=False, type=float, default=0.05, help='Maximum minor allele frequency value to filter variants')
    parser.add_argument('-t', '--threads', required=False, type=int, default=48, help='Threads to run in parallel')
    return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
