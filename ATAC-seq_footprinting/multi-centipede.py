#!/usr/bin/python3

import os
import sys
import subprocess
from multiprocessing import Pool

sample = sys.argv[1]

def runCentipede(motifName, sampleName):
    d_base = '/home/joshchiou/joshchiou-data2/FACS.alpha+beta_cells/'
    d_bam = os.path.join(d_base, 'bams')
    d_bed = os.path.join(d_base, 'motifs', sampleName)
    d_matrix = os.path.join(d_base, 'matrices', sampleName)
    d_centipede = os.path.join(d_base, 'centipede', sampleName)
    d_log = os.path.join(d_centipede, 'log')
    d_post_centipede = os.path.join(d_base, 'post-centipede', sampleName)   

    readsBam = os.path.join(d_bam, '.'.join([sampleName, 'bam']))
    motifBed = os.path.join(d_bed, '.'.join([motifName, sampleName, 'bed']))
    motifMatrix = os.path.join(d_matrix, '.'.join([motifName, sampleName, 'discrete.matrix.gz']))
    centipedeLog = os.path.join(d_centipede, 'log', '.'.join([motifName, sampleName, 'centipede.log'])) 
    
    d_list = [d_base, d_bam, d_bed, d_matrix, d_centipede, d_log, d_post_centipede]
    for d in d_list:
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass
    if os.path.exists(motifBed) and os.path.getsize(motifBed) > 0:
        try:
            with open(centipedeLog, 'w') as f:
                subprocess.call(
                        'make_cut_matrix -d -b \'(1-100000 1)\' -p 1 {0} {1} | gzip -c > {2}'.format(readsBam, motifBed, motifMatrix), 
                        stdout=f,
                        stderr=f,
                        shell=True)
        except:
            sys.stderr.write(Exception)
            sys.stderr.write('\n')
            pass

    if os.path.isfile(motifMatrix):
        try:
            with open(centipedeLog, 'a') as f:
                subprocess.call('./centipede.R {0} {1} {2}'.format(motifMatrix, motifBed, sampleName), 
                        stdout=f,
                        stderr=f,
                        shell=True)
        except:
            sys.stderr.write(Exception)
            sys.stderr.write('\n')
            pass
    
    centipedeOut = os.path.join(d_centipede, '.'.join([motifName, sampleName, 'PostPr.txt']))
    combinedOut = os.path.join(d_post_centipede, '.'.join([motifName, sampleName, 'footprints.bed']))
    
    if os.path.isfile(centipedeOut):
        subprocess.call('paste {0} {1} | awk \'$7 > 0.95\' > {2}'.format(motifBed, centipedeOut, combinedOut), shell=True)
    return

with open('/home/joshchiou/references/combined_motifs.list','r') as f:
    motifList = f.read().splitlines()

argsList=[]
for motif in motifList:
    argsList.append((motif, sample))

pool = Pool(processes=(24))
pool.starmap(runCentipede, argsList)
pool.close()
pool.join()
