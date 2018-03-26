.rs.restartR()
library(devtools)
library(GenomicRanges)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(Homo.sapiens)
library(doMC)
library(plyr)
registerDoMC(cores = 16) # The number of cpu-cores used in the analysis
options(warn=-1)

# devtools::install_github(repo = "joe-nas/haploplotR")

devtools::build()
devtools::load_all()

data(analysis1dat)

ldabase <- LDABase$new(file_path = "/home/SSD-Data/1000Genomes/")
ldabase <- LDABase$new(file_path = "/scratch/jfalck/1000G/")

# This sets up an object specifying the analysis interval and which populations are to be analyzed
res2 <- llply(analysis_intervals_gr[1:10], function(x){
  dat <- LDAImport$new(ldabase = ldabase, granges =  x, populations = populations)$set_data()
  res <- LDAanalysis1$new(lead_snps = lead_snps, cutoff = 0.8, lda_import = dat)
  res$set_rsquared()
  res$set_results()
  invisible(res)
}, .parallel = F, .inform = F, .progress = "text", warn.conflicts = FALSE)


# this collects tag snp intervals from multiple regions into one object
high_ld2 <- Reduce(c,Map(function(x) x$results, res2))
save(high_ld2, file = "data/high_ld2_analysis1.RData")
