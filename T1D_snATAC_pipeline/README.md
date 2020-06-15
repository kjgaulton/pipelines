# Instructions
Run the `snATAC_pipeline_CB.py` (for combinatorial indexing single cell ATAC-seq)  or `snATAC_pipeline_10X.py` (for droplet-based 10X single cell ATAC-seq) pipeline to generate input matrices for clustering.

Use the `T1D_merged_snATAC.ipynb` notebook using the input matrices to cluster cells.

# Requirements
## Software packages  
R    
trim_galore  
bwa  
samtools  
bedtools  
picard  
macs2  

## python3 packages
numpy  
pandas  
scipy  
matplotlib  
seaborn  
pysam  
rpy2  
scanpy  
anndata  

# R packages
Matrix  
harmony  
