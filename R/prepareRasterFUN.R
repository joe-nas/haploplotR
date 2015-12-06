#' Maps numeric values of matrix to colors
#'
#' @param mat E.g. R^2 Symmetric numeric matrix with range [0,1] representing pairwise R^2 values.
#' @param depth The  depth of pairwise LD to be shown.
#' @param colorscheme
#' @param mapping DataFrame with columns 'pos' and 'ids' specifying SNP positions
#' and SNP names.
#' @param color_highlights List specifying highlight colors in mapping panel. Must have
#' common, lead and tag slots.
#' @return data.frame The same as mapping but with tag and alpha columns
#' @examples
#' color_highlights <- list(common = 'green', lead = 'yellow', tag = 'cyan')

prepareRaster <- function(mat, depths, colorscheme){
  dms <- nrow(mat)
  mask <- upper.tri(mat, diag = F)
  if(depths > dms/2){
    mask[1:depths, depths:dms] <- lower.tri(mask[1:depths, depths:dms],diag = F) #depths > dms/2
  }else{
    mask[1:(dms-depths), depths:dms] <- lower.tri(mask[1:(dms-depths), depths:dms],diag = F) #depths < dms/2
    mask[lower.tri(mask, diag = T)] <- FALSE
  }

  mat[mask == T] <- colorscheme[mat[mask == T] * 100 + 1]
  mat[mask == F] <- "transparent"
  mat
}
