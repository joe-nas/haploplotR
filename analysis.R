library(devtools)
library(GenomicRanges)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
devtools::load_all('/home/SSD-Data/Projects/haploplotR/')

lead_snps <- c("rs12569773", "rs4615961", "rs10983654", "rs1512262", "rs11119348" )
populations <- c("CEU", "CHB", 'ABC')


