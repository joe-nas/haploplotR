library(R6)

LDABase <- R6Class("LDABase",
  public = list(
    chromosome = NULL,
    start = NULL,
    end = NULL,
    position = NULL,
    lead_snps = NULL,
    initialize = function(chromosome = NA, start = NA, end = NA, lead_snps = NA){
      self$chromosome <- chromosome
      self$start <- start
      self$end <- end
      self$position <- GRanges(Rle(chromosome), IRanges(start = start, end = end))
      self$lead_snps <- lead_snps
    }
  )
)

LDAPopulation <- R6Class("LDApopulation",
                         inherit = LDABase,
                         public = list(
                           population = NULL,
                           data = NULL,
                           info = NULL,
                           rsquared = NULL,
                           vcf_file = NULL,
                           panel_file = NULL,
                           initialize = function(population = NA, vcf_file = "",
                                                 panel_file = ""){
                             self$population  <- population
                           }
                         ))


LDAnalysis <- R6Class("LDAnalysis",
                  public = list(
                    chromosome = NULL,
                    startp = NULL,
                    endp = NULL,
                    position = NULL,
                    lead_snps = NULL,
                    population = NULL,
                    data = NULL,
                    info = NULL,
                    vcf_file = NULL,
                    panel_file = NULL,
                    initialize = function(chromosome = NA, startp = NA, endp = NA, lead_snps = NA, population = NA, vcf_file = NA,
                                          panel_file = NA) {
                      self$chromosome <- chromosome
                      self$startp <- startp
                      self$endp <- endp
                      #self$position <- GRanges(Rle(chromosome), IRanges(start = startp, end = endp))
                      self$lead_snps <- lead_snps
                      self$population <- population
                      self$vcf_file <- vcf_file
                      self$panel_file <- panel_file
                      #self$greet()
                    },
                    set_postition = function() {
                      self$position <- GRanges(Rle(self$chromosome), IRanges(start = self$startp, end = self$endp))
                    },
                    set_data = function(where = self$position, which = self$population, vcf_file = self$vcf_file,
                                        panel_file = self$panel_file) {
                      self$data <- get_SnpMatrix_from_vcf(where = where, which = which, vcf_file = vcf_file,
                                                          panel_file = panel_file)
                    },
                    greet = function() {
                        cat(paste0("Analysis takes place in ", self$position, ".\n"))
                    }
                  )
)

https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html

lda1 <- LDAnalysis$new(chromosome = "chr1" , startp = 122831384, endp = 123481383, population = "CEU",panel_file = "integrated_call_samples_v3.20130502.ALL.panel", vcf_file = "/home/SSD-Data/1000Genomes/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

