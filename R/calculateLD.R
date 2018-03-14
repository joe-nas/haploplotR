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
  idx <- which(colnames(xSnpMatrix) %in% lead_snps)
  if(length(idx) != 0){
    ySnpMatrix <- xSnpMatrix[,which(colnames(xSnpMatrix) %in% lead_snps)] # limit the calculation to this lead_snps
    snpStats::ld(x = xSnpMatrix, y = ySnpMatrix, stats = c('R.squared'))
  }else{
    NULL
  }

  # this function needs optimization! when lead_snps are not found in the data, pairwise LD is calculated on xSnpMatrix x xSnpMatrix
#  if(dim(ySnpMatrix == xSnpMatrix)){
#    colnames(xSnpMatrix)
#  }else{
#  }
}

## debug
#calculateLD(lead_snps,TGData[[11]]$data$LWK$genotype)
