#!/usr/bin/python3

import os
import sys
import subprocess
from multiprocessing import Pool

sample = sys.argv[1]

def runFimo(motifName, sampleName):
	d_base = '/home/joshchiou/joshchiou-data2/FACS.alpha+beta_cells/'
	d_seqs = os.path.join(d_base, 'seqs')
	d_bed = os.path.join(d_base, 'motifs', sampleName)

	dataPWM = os.path.join('/home/joshchiou/references/combined.motif_database.meme')
	seqFasta = os.path.join(d_seqs, '.'.join([sampleName, 'seqs.fa']))
	motifBed = os.path.join(d_bed, '.'.join([motifName, sampleName, 'bed']))

	d_list = [d_base, d_seqs, d_bed]
	for d in d_list:
		if not os.path.exists(d):
			try:
				os.makedirs(d)
			except:
				pass

	try:
		with open(os.devnull, 'w') as f:
			subprocess.call('fimo --parse-genomic-coord --max-strand --skip-matched-sequence --bgfile motif-file --motif {0} {1} {2} | awk \'BEGIN{{FS=OFS=\"\\t\"}} NR>1 {{print $3,$4,$5+1,$1,$7,$6}}\' | sort -k1,1 -k2,2n | awk -F\"\\t\" \'!uniq[$1 FS $2 FS $3 FS $6]++\' > {3}'.format(motifName, dataPWM, seqFasta, motifBed), stderr=f, shell=True)
	except:
		sys.stderr.write(Exception)
		sys.stderr.write('\n')
		pass
	return

with open('/home/joshchiou/references/combined_motifs.list','r') as f:
	motifList = f.read().splitlines()

argsList=[]
for motif in motifList:
	argsList.append((motif, sample))

pool = Pool(processes=(24))
pool.starmap(runFimo, argsList)
pool.close()
pool.join()
