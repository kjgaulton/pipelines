#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import subprocess
import logging
from multiprocessing import Pool


def process_vcf(args, footprints, name):
    logging.info('Processing motif: {}'.format(name))
    vcf_output = os.path.join(args.output, name + '.vcf')
    with open(footprints) as f, open(vcf_output, 'w') as vcf_out:
        lines = f.read().splitlines()
        tmp = tempfile.TemporaryFile()
        with open(args.vcf_header) as head:
            vcf_out.write(head.read())
        for line in lines:
            chrom, start, end, motif = line.split('\t')
            vcf_file = os.path.join(args.vcf_directory, 'GoT2D.chr{}.paper_integrated_snps_indels_sv_beagle_thunder.vcf.gz'.format(chrom)) 
            query = '{0}:{1}-{2}'.format(chrom, start, end)
            tabix_cmd = ['tabix', '-f', vcf_file, query]
            subprocess.call(tabix_cmd, stdout=tmp)
        tmp.seek(0)
        records = [rec.decode() for rec in tmp.read().splitlines()]
        for rec in records:
            if len(rec.split('\t')[3]) > 1: continue
            if len(rec.split('\t')[4]) > 1: continue
            if rec.split('\t')[6] != 'PASS': continue
            info = rec.split('\t')[7]
            AC = float(info.split(';')[0].split('=')[1])
            AN = float(info.split(';')[1].split('=')[1])
            RSQ = float(info.split(';')[6].split('=')[1])
            MAF = AC/AN
            if MAF > args.maximum_maf and 1-MAF > args.maximum_maf: continue
            if RSQ < args.minimum_rsquared: continue
            vcf_out.write('{}\n'.format(rec))
        tmp.close()
    bgzip_cmd = ['bgzip', vcf_output]
    index_cmd = ['tabix', vcf_output + '.gz']
    if os.path.exists(vcf_output) and os.path.getsize(vcf_output) != 0:
        subprocess.call(bgzip_cmd)
    if os.path.exists(vcf_output+'.gz'):
        subprocess.call(index_cmd)
    return

def main(args):
    logging.info('Starting up.')
    with open('fp.list') as f:
        fp_names = f.read().splitlines()
    fp_files = ['merged-footprints/' + name + '.bed' for name in fp_names]
    args_list = list(zip([args for i in fp_names], fp_files, fp_names))
    pool = Pool(processes=(args.threads))
    pool.starmap(process_vcf, args_list)
    pool.close()
    pool.join()
    logging.info('Finishing up.')
    return

def process_args():
    parser = argparse.ArgumentParser(description='Make VCF files for each footprint motif.')
    #parser.add_argument('-i', '--footprints', required=True, type=str, help='Path to footprints file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to output directory')
    #parser.add_argument('-n', '--name', required=True, type=str, help='Output name to prepend to file')
    parser.add_argument('-vcf_dir', '--vcf_directory', required=True, type=str, help='Path to directory containing tabix-indexed vcfs')
    parser.add_argument('-head', '--vcf_header', required=False, type=str, default='etc/header.txt', help='Path to vcf header')
    parser.add_argument('-min_rsq', '--minimum_rsquared', required=False, type=float, default=0.6, help='Minimum R^2 value to filter variants')
    parser.add_argument('-max_maf', '--maximum_maf', required=False, type=float, default=0.05, help='Maximum minor allele frequency value to filter variants')
    parser.add_argument('-t', '--threads', required=False, type=int, default=48, help='Threads to run in parallel')
    return parser.parse_args()

args = process_args()
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
main(args)
