




identifyHighLD <- function(rsquared, info, cutoff = 0.8, population = NULL){
  #identify tag-SNPs
  res <- which(rsquared >= cutoff, arr.ind = T) # rownames: all SNPs in ranges, colnames lead-SNPS
  lead_snps_v <- colnames(rsquared)[res[,2]] # lead_snp vector
  res_info_df <- data.frame(lead_snp = lead_snps_v, info[res[,1],], rsquared = rsquared[res])
  chromosomes <- Rle(paste("chr",res_info_df$chromosome_num, sep=""))
  tag_snp_r <- as.numeric(rownames(res_info_df)) #tag sno rows
  startP <-info[tag_snp_r-1,]$pos # start position
  endP <- info[tag_snp_r+1,]$pos # end position
  # results
  tag_snps <- GRanges(chromosomes, IRanges(star = startP, end = endP),
                      mcols = data.frame(lead_pos = res_info_df$pos, lead_snp = res_info_df$lead_snp,
                                         tag_snp = res_info_df$ids, rs = res_info_df$rsquared, population = population))
  tag_snps
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



