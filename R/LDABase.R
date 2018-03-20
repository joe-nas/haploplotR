require(R6)
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
