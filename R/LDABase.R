#' Object for holding and file needed for the calculation of linkage disequilibrium
#'
#' @param lead_snps vector contain dbsnpids
#' @param populations vecotr containing 1000genomes population abbreviations
#' @param file_path working directory
#' @param vcf_file vcf filename
#' @param panel_file panel filename
#'
#' @return
#' @export
#'
#' @examples

require(R6, quietly = T, warn.conflicts = FALSE)
LDABase <- R6Class("LDABase",
                   public = list(
                     lead_snps = NULL,
                     populations = NULL,
                     vcf_file = NULL,
                     panel_file = NULL,
                     file_path = NULL,
                     initialize = function(lead_snps = NA, populations = NA, file_path = ".",
                                           vcf_file = "ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
                                           panel_file = "integrated_call_samples_v3.20130502.ALL.panel"){
                       self$lead_snps <- lead_snps
                       self$populations <- populations
                       self$panel_file <- file.path(file_path,panel_file)
                       self$vcf_file <- file.path(file_path, vcf_file)
                       self$file_path <- file_path
                     }
                   ))
