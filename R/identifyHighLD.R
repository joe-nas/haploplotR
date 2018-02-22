




identifyHighLD <- function(rsquared, info, cutoff = 0.8){
  res <- which(rsquared >= cutoff, arr.ind = T) # rownames: all SNPs in ranges, colnames lead-SNPS
  res_info <- data.frame(info[res[,1],],rsquared = rsquared[res])
  res_info
}



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



