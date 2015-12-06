#' Preparing mapping data.frame in upperTriGrob
#'
#' @param r2 Symmetric numeric matrix with range [0,1] representing pairwise R^2 values.
#' @param lead Vector of SNP identifiers for which tag SNPs are identifierd and
#' which are be highlighted in the mapping panel.
#' @param r2_cutoff The cutoff value used to identify tag SNPs.
#' @param mapping DataFrame with columns 'pos' and 'ids' specifying SNP positions
#' and SNP names.
#' @param color_highlights List specifying highlight colors in mapping panel. Must have
#' common, lead and tag slots.
#' @return data.frame The same as mapping but with tag and alpha columns
#' @examples
#' color_highlights <- list(common = 'green', lead = 'yellow', tag = 'cyan')


prepareMapping <-
  function(r2, lead, r2_cutoff, mapping, color_highlights = list(
    common = 'black', lead = 'blue', tag = 'red')){

  lead <- unique(lead[lead %in% rownames(r2)])
  mapping$color <- color_highlights$common
  mapping$alpha <- 0.2
  if(!length(lead) == 0){
    tag <- which(r2[lead,] >= r2_cutoff, arr.ind = T)
    tag <- if(length(lead) > 1){
      tag[,2]
    }else{
      tag
    }
    mapping[tag,'color'] <- color_highlights$tag
    mapping[mapping$ids %in% lead,'color'] <- color_highlights$lead
    mapping[tag,'alpha'] <- 1
  }
  invisible(mapping)
}
