#' From SnpMatrix objects generated with import1000GData function generade LD data
#'
#' @param
#' @param lead_snps the lead_snps to be analyzed
#' @param vcf_file local file or url to vcf file.
#' @param panel_file local file or url to panel file.
#' @return Object of class SnpMatrix.
#' @examples
#' COMPLETE EXAMPLES




calculateLD <- function(lead_snps, xSnpMatrix, ...){
  ySnpMatrix <- xSnpMatrix[,which(colnames(xSnpMatrix) %in% lead_snps)] # limit the calculation to this lead_snps

  # this function needs optimization! when lead_snps are not found in the data, pairwise LD is calculated on xSnpMatrix x xSnpMatrix
#  if(dim(ySnpMatrix == xSnpMatrix)){
#    colnames(xSnpMatrix)
#  }else{

    snpStats::ld(x = xSnpMatrix, y = ySnpMatrix, stats = c('R.squared'))

#  }
}

my_mapping <- Reduce(rbind,lapply(seq_along(target_dat), function(i){
  prepareMapping(r2 = target_rsquared[[i]],mapping =  target_dat[[i]]$info,r2_cutoff = 0.8,
                 lead = c('rs10863790','rs10863790','rs11119348', 'rs12566152','rs861020','rs10863790'))
}))
