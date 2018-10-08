# ChIP-seq imbalance pipeline

This vignette covers the three steps of ChIP-seq imbalance analysis (on
gatsby):

1. genotyping (QuASAR)
1. remapping (WASP)
1. computing imbalance statistics

We will apply these steps to two groups of FASTQ files: islet sample HI-87 from
the pasquali data, and our recent 2018 HepG2/FOXA2 data.

| Dataset | No. FASTQ files | Details |
|---|---|---|
| HI-87 | 2 | Chip-seq of MAFB and NKX2-2 in islet |
| HepG2 FOXA2 | 3 | 3x ChIP-seq of FOXA2 in HepG2 |

The first section (Example Runs) shows how these datasets are processed in a
dense format. The subsequent sections break down each step in more detail.

## Example Runs

### HI-87

We can process the HI_87 data (two datasets for two different transcription
factors) with these commands:

```sh
cd /home/data/chip-imbalance-example/
quasar-genotype raw-reads/*_HI_87.fastq.gz \
  --write-bam \
  --bam-dir alignment \
  --vcf-chr genotype/HI_87 \
  --sample HI_87 \
  --quiet \
  --processes 8
for vcf in genotype/HI_87.chr*.vcf; do bgzip $vcf; done
wasp-map remap \
  --sample MAFB_HI_87 alignment/ERS353636_MAFB_HI_87.filt.bam \
  --sample NKX2_2_HI_87 alignment/ERS353636_NKX2_2_HI_87.filt.bam \
  --vcf-format genotype/HI_87.chr{}.vcf.gz \
  --vcf-sample HI_87 \
  --r2 0 \
  --processes 16
for factor in MAFB NKX2_2
do pileup-to-counts --json remap/pileup/${factor}_HI_87.pileup \
  | call-imbalance --cores 16 \
  > ${factor}_HI_87_imbalance.txt
done 
```

### HepG2 FOXA2

We can process the 2018 HepG2 FOXA2 data (three datasets for the same cell type
and transcription factor) with these commands:

```sh
cd /home/data/chip-imbalance-example/
quasar-genotype raw-reads/HepG2_FOXA2_*_L004_R1_001.fastq.gz \
  --write-bam \
  --bam-dir alignment \
  --vcf-chr genotype/HepG2_FOXA2 \
  --sample HepG2_FOXA2 \
  --quiet \
  --processes 8
for vcf in genotype/HepG2_FOXA2.chr*.vcf; do bgzip $vcf; done
wasp-map remap \
  --sample HepG2_FOXA2_1 alignment/HepG2_FOXA2_1_S1_L004_R1_001.filt.bam \
  --sample HepG2_FOXA2_2 alignment/HepG2_FOXA2_2_S2_L004_R1_001.filt.bam \
  --sample HepG2_FOXA2_3 alignment/HepG2_FOXA2_KG02_Human_M_23_S2_L004_R1_001.filt.bam \
  --vcf-format genotype/HepG2_FOXA2.chr{}.vcf.gz \
  --vcf-sample HepG2_FOXA2 \
  --r2 0 \
  --processes 16
cluster-pileups --processes 16 remap/pileup/HepG2_FOXA2_*.pileup \
  | call-imbalance --cores 16 \
  > HepG2_FOXA2_imbalance.txt
```

## Step 1: Genotyping

The genotyping step is performed by a python script at
`/home/data/aaylward-source/quasar-genotype/quasar_genotype.py`. It is symlinked
to `/usr/local/bin/quasar-genotype`, so it can be called simply with
`quasar-genotype`. For example, to get its help text:

```sh
quasar-genotype -h | less
```

There are a lot of options, but most have good defaults so in practice only
a few need to be considered. The simplest call requires only the sequencing
data files, which can be in FASTA, FASTQ, or BAM format. Here is an example
using islet sample HI_87:

```sh
cd /home/data/chip-imbalance-example/
quasar-genotype raw-reads/ERS353636_MAFB_HI_87.fastq.gz raw-reads/ERS353636_NKX2_2_HI_87.fastq.gz
```

If there are more than a few data files, a glob pattern can be used to make
things easier

```sh
quasar-genotype raw-reads/*_HI_87.fastq.gz
```

This will align the reads, infer genotypes, and print a VCF file to standard
output. That doesn't provide exactly the inputs required by the next step
(remapping with WASP) though, so a more appropriate command is:

```sh
quasar-genotype raw-reads/*_HI_87.fastq.gz \
  --write-bam \
  --bam-dir alignment \
  --vcf-chr genotype/HI_87 \
  --sample HI_87 \
  --quiet \
  --processes 8
```

As a follow-up, the resulting VCF files need to be compressed before submitting
them for imputation or using them with WASP.

```sh
for vcf in genotype/HI_87.chr*.vcf; do bgzip $vcf; done
```

### Parameter Details

The `--write-bam` option means the aligned reads will be written to disk in BAM
format. The `--bam-dir` option indicates that the BAM files will be placed in
the `alignment` directory, and the `--vcf-chr` option specifies a prefix for
writing VCF files split by chromosome. The `--sample` option sets the sample
name in the VCF header, in this case to `HI_87`, and the `--quiet` option
prevents writing a VCF to standard output.

### Results

After running this command, we can see the BAM files written to the `alignment`
directory:

```sh
ls alignment/*HI_87*
```

```
ERS353636_MAFB_HI_87.align.log  ERS353636_MAFB_HI_87.filt.bam.bai  ERS353636_NKX2_2_HI_87.filt.bam
ERS353636_MAFB_HI_87.filt.bam   ERS353636_NKX2_2_HI_87.align.log   ERS353636_NKX2_2_HI_87.filt.bam.bai
```

and the VCF files written to the `genotype` directory:

```sh
ls genotype/HI_87.chr*.vcf
```

```
HI_87.chr10.vcf  HI_87.chr13.vcf  HI_87.chr16.vcf  HI_87.chr19.vcf  HI_87.chr21.vcf  HI_87.chr3.vcf  HI_87.chr6.vcf  HI_87.chr9.vcf
HI_87.chr11.vcf  HI_87.chr14.vcf  HI_87.chr17.vcf  HI_87.chr1.vcf   HI_87.chr22.vcf  HI_87.chr4.vcf  HI_87.chr7.vcf  HI_87.chrX.vcf
HI_87.chr12.vcf  HI_87.chr15.vcf  HI_87.chr18.vcf  HI_87.chr20.vcf  HI_87.chr2.vcf   HI_87.chr5.vcf  HI_87.chr8.vcf
```

```sh
zcat genotype/HI_87.chr1.vcf.gz | head
```

```
##fileformat=VCFv4.0
##reference=GRCh37
##INFO=<ID=BUILD,Number=1,Type=Integer,Description="Genome build">
##INFO=<ID=GP,Number=3,Type=Float,Description="Genotype probabilities">
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HI_87	HI_87_COMP
1	740416	rs116024440	A	G	.	PASS	BUILD=37;GP=0.999874454972485,0.000125545027514821,4.59405542148581e-35	GT	0/0	0/1
1	787399	rs2905055	G	T	.	PASS	BUILD=37;GP=2.94476993103003e-34,0.00101032607630337,0.998989673923697	GT	1/1	0/0
1	794172	rs11804026	G	T	.	PASS	BUILD=37;GP=0.997642876228381,0.00235712377161868,3.84716962751783e-22	GT	0/0	0/1
1	794332	rs12127425	G	A	.	PASS	BUILD=37;GP=0.999939168135806,6.08318641942384e-05,6.85367077069993e-41	GT	0/0	0/1
```

## Step 2: Remapping

The WASP pipeline script is at `/home/data/aaylward-source/wasp-map/wasp_map.py`
It is symlinked to `/usr/local/bin/wasp-map` so it can be called simply with
`wasp-map`. For example:

```sh
wasp-map -h | less
```

To perform the remapping, we need only provide the BAM and VCF files from
`quasar-genotype` as inputs to `wasp-map`:

```sh
wasp-map remap \
  --sample MAFB_HI_87 alignment/ERS353636_MAFB_HI_87.filt.bam \
  --sample NKX2_2_HI_87 alignment/ERS353636_NKX2_2_HI_87.filt.bam \
  --vcf-format genotype/HI_87.chr{}.vcf.gz \
  --vcf-sample HI_87 \
  --r2 0 \
  --processes 16
```

### Parameter Details

The first argument is an output directory, in this case `remap`. The `--sample`
options indicate a sample name and a file path for each BAM input file. The
`--vcf-format` option indicates the form of the input VCF files (The braces
`{}` will be replaced with a chromosome number from 1 to 22). The
`--vcf-sample` option indicates the sample ID on the input VCF files (in this
case `HI_87`, just as specified in step 1). Here, the `--r2` option is used
to set the imputation quality threshold to 0, since imputed genotype data is
not being used. If the input VCFs are from Michigan Imputation Server, the
`--r2` option can be omitted and the pipeline will use a default threshold of
0.9.

### Results

After running this command, the results of interest (pileups) can be found in
the `pileup` subdirectory of the output directory.

```sh
head remap/pileup/*_HI_87.pileup
```

```
==> remap/pileup/MAFB_HI_87.pileup <==
chr1	1098421	C	0	*	*
chr1	1310924	T	1	c	J
chr1	1813237	C	2	,,	>J
chr1	2478511	T	1	a	B
chr1	7158287	C	4	t.T.	HJEF
chr1	7322807	A	3	G.G	JIJ
chr1	7684030	G	1	A	<
chr1	8320045	G	8	.a.aa,,a	DHJJGJGJ
chr1	8430594	C	3	.$..	8JI
chr1	8431105	C	2	.,	IJ

==> remap/pileup/NKX2_2_HI_87.pileup <==
chr1	1098421	C	15	,$t,T,t,t,,ttt,,	?1F.F?JJJGHFGJE
chr1	1310924	T	12	.,,,,Cccc,c,	IHHJJJJJJIGJ
chr1	1342612	G	6	.CCC.C	JGIGFJ
chr1	1407043	A	9	,CcC.c,c,	FDHFHDHDD
chr1	1813237	C	6	t,ttT,	HIJJJJ
chr1	2478511	T	5	,,,aa	HDJIJ
chr1	2574078	T	6	C,.c,.	JHJHFF
chr1	4587880	G	8	,a.,.a,.	JJGIJGA@
chr1	5373029	t	6	.ccC..	IFHJFC
chr1	6614535	G	6	A...AA	1?@D<9
```

## Step 3: Computing imbalance statistics

This step involves three scripts: `pileup-to-counts`, `cluster-pileups`, and
`call-imbalance`. Like the previous scripts, each has a help message viewable
via the `-h` flag.

```sh
pileup-to-counts -h 
cluster-pileups -h
call-imbalance -h
```

Also like the previous scripts, these are symbolic links to the scripts' actual
locations:

```
/home/data/aaylward-source/pileup-to-counts/pileup_to_counts.py
/home/data/aaylward-source/cluster-pileups/cluster_pileups.py
/home/data/aaylward-source/call-imbalance/call_imbalance.R
```

How they are applied to the HI-87 data differs from how they are applied to the
HepG2 FOXA2 data.

### HI-87

With respect to the HI-87 data, we are ready to compute the imbalance
statistics using `pileup-to-counts` and `call-imbalance`:

```sh
for factor in MAFB NKX2_2
do pileup-to-counts --json remap/pileup/${factor}_HI_87.pileup \
  | call-imbalance --cores 16 \
  > ${factor}_HI_87_imbalance.txt
done 
```

`pileup-to-counts` converts a pileup into allele counts that can be processed
by `call-imbalance`.

Here the `--json` option causes `pileup-to-counts` to produce the allele count
data as a JSON file instead of a normal tabular file. This is because
`call-imbalance` is an R script and I've found that R isn't great at reading
tabular data through a pipe.

### HepG2 FOXA2

In the case of the HepG2 FOXA2 data, this step differs.
Rather than one dataset per transcription factor, we have a handful of FOXA2
datasets that we wish to combine and analyze together. For this purpose, we can
use `cluster-pileups`:

```sh
cluster-pileups --processes 16 remap/pileup/HepG2_FOXA2_*.pileup \
  | call-imbalance --cores 16 \
  > HepG2_FOXA2_imbalance.txt
```

`cluster-pileups` estimates an empirical beta-binomial distribution for each
input pileup, identifies clusters of pileups with similar distributions,
merges each cluster into one set of count data, and then passes these count
datasets to `call-imbalance` to handle statisticslly.

### Results

Each output file contains some statistics quantifying allelic imbalance at
each variant that passes a minimum coverage threshold (10-fold, by default).

```sh
head MAFB_HI_87_imbalance.txt 
```
```
chr     pos     coverage        ref_count       binom_pval      binom_fdr       beta_binom_pval beta_binom_fdr  log_posterior_allelic_fold_change
1       163064483       16      9       1       1       1       1       0.150808631087479
1       163064579       15      6       0.300461912299701       1       0.405738282511011       1       -0.44089224755327
1       163064766       13      6       0.582399317614533       1       0.720599498089526       1       -0.205848755427337
1       208698573       11      5       0.555506684067047       1       0.737521773646807       1       -0.213806639923069
1       227063400       12      7       1       1       1       1       0.204019376880525
3       2519703 11      8       0.365320859141406       1       0.440303683264956       1       0.671571180666862
3       111732608       10      6       1       1       1       1       0.239146614737931
3       121865446       10      4       0.357676938376882       1       0.541862757563377       1       -0.370829428120181
3       133132496       11      5       0.555506684067047       1       0.737521773646807       1       -0.213806639923069
```

I'm calling the effect size column "log posterior allelic fold change," since it's similar to log allelic fold change but is based on a posterior estimate of the "ground truth" allelic fraction (computed from the observed counts and the estimated beta prior), rather than just the observed counts.
