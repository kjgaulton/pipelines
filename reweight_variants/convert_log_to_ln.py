#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

metabo = pd.read_csv(sys.argv[1], sep='\t', header=None)
metabo.columns = ['chrom', 'position', 'rsid', 'locus', 'logbf', 'state']
base10 = np.array([10 for x in metabo['logbf']])
metabo['lnbf'] = np.log(np.power(base10, metabo['logbf']))

out_cols = ['chrom', 'position', 'rsid', 'locus', 'lnbf', 'state']
metabo_out = metabo[out_cols]
metabo_out.to_csv('metabo.islet_ATAC.lnbf.bed', sep='\t', header=False, index=False)
