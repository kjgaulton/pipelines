#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import logging
from multiprocessing import Pool


def make_groupfiles(args, name):
    try:
        vcf_file = os.path.join(args.proc_vcf_dir, name + '.vcf.gz')
        footprints_file = os.path.join(args.footprints_dir, name + '.bed')
        p1_cmd = 'zcat {}'.format(vcf_file)
        p2_cmd = """awk '!/^#/'"""
        p3_cmd = 'cut -f -5'
        p4_cmd = """awk 'BEGIN{FS=OFS="\t"} {print $1,$2-$1,$2,$1":"$2"_"$4"/"$5}'"""
        p5_cmd = 'bedtools intersect -a - -b {} -wa -wb'.format(footprints_file)
        p_cmd = '{} | {} | {} | {} | {}'.format(p1_cmd, p2_cmd, p3_cmd, p4_cmd, p5_cmd)
        intersect_out = os.path.join(args.intersect_dir, name + '.var-fp.bed')
        with open(intersect_out, 'w') as f:
            subprocess.call(p_cmd, stdout=f, shell=True)
        
        disrupt_out = os.path.join(args.disrupt_dir, name + '.disrupt.out')
        disrupt_err = os.path.join(args.disrupt_dir, name + '.disrupt.err')
        with open(disrupt_out, 'w') as out, open(disrupt_err, 'w') as err:
            p6 = subprocess.Popen(['python3', 'motif_disruption.py', intersect_out], stdout=subprocess.PIPE, stderr=err)
            p7 =subprocess.call(['awk', r'$14 <= 1.0'], stdin=p6.stdout, stdout=out)
        
        group_file = os.path.join(args.group_dir, name + '.group')
        with open(group_file, 'w') as f, open(disrupt_out) as out:
           variants = [line.split('\t')[3] for line in out.read().splitlines()]
           variants = sorted(list(set(variants)))
           group_line = [name] + variants
           f.write('\t'.join(group_line))
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
    parser = argparse.ArgumentParser(description='Make group files for EPACTS run.')
    parser.add_argument('--proc_vcf_dir', required=True, type=str, help='Path to processed VCFs directory')
    parser.add_argument('--footprints_dir', required=True, type=str, help='Path to unmerged footprints directory')
    parser.add_argument('--intersect_dir', required=True, type=str, help='Path to variants-footprints intersect directory')
    parser.add_argument('--disrupt_dir', required=True, type=str, help='Path to disrupted variants directory')
    parser.add_argument('--group_dir', required=True, type=str, help='Path to group files directory')
    parser.add_argument('-t', '--threads', required=False, type=int, default=48, help='Threads to run in parallel')
    return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
