




identifyHighLD <- function(rsquared, info, cutoff = 0.8){
  #identify tag-SNPs
  res <- which(rsquared >= cutoff, arr.ind = T) # rownames: all SNPs in ranges, colnames lead-SNPS
  res_info_df <- data.frame(info[res[,1],],rsquared = rsquared[res])
  chromosomes <- Rle(paste("chr",res_info_df$chromosome_num, sep=""))
  startP <-info[as.numeric(rownames(res_info_df))-1,]$pos
  endP <- info[as.numeric(rownames(res_info_df))+1,]$pos
  tag_snps <- GRanges(chromosomes,IRanges(star = startP, end = endP), mcols = res_info_df)
  tag_snps
}

constructTagSNPIntervals <-function(res_info){
  tag_SNP_rowpos <- as.integer(rownames(res_info))

}

  #identifyHighLD(TGData[[3]]$rsquared$CEU, TGData[[3]]$data$CEU$info )

# tagSNPintervals <- function(){
#
# }
#
#
# results_df <- data.frame(chromosome=vector(),start=vector(),end=vector(),tag_snp=vector(),lead_snp=vector(),rsquared=vector())
#
# generate
#
#
# analysis1object
#   - granges
#     - list(populations)
#       - data (1000G)



