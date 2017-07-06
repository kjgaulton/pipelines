#!/usr/bin/Rscript
# centipede.R
# Usage: /usr/bin/Rscript --vanilla centipede.R matrix.file motifs.file

#dir.base="/data2/smorabito-data2/recreate-ATAC-footprints"

# start of script
library(CENTIPEDE)
library(tools)

# process args
matrix.in <- commandArgs()[7]
motifs.in <- commandArgs()[8]
print(paste0("Input matrix: ", matrix.in))
print(paste0("Input motifs: ", motifs.in))
file.name <- basename(file_path_sans_ext(motifs.in))

# output centipede posterior probabilities
#dir.centipede <- file.path(dir.base,"centipede")
#dir.output.centipede <- file.path(dir.base, file.name)
#setwd(dir.output.centipede)

num.columns <- max(count.fields(gzfile(matrix.in), sep="\t"))
count.matrix <- read.table(gzfile(matrix.in), sep="\t", header=F, fill=T, col.names=paste0('V', seq_len(num.columns)))
count.matrix[is.na(count.matrix)] <- 0
motifs <- read.table(motifs.in, sep="\t", header=F)
centFit <- fitCentipede(Xlist = list(ATAC=as.matrix(count.matrix)), Y=cbind(rep(1, dim(motifs)[1]), motifs[,5]))
write.table(centFit$PostPr,col=F,row=F,sep="\t",quote=F,file=paste0(basename(matrix.in),".txt"))
#write(centFit$PostPr, stdout())