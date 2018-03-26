require(R6, quietly = T)
LDAImport <- R6Class("LDAImport",
                     public = list(
                       ldabase = NULL,
                       granges = NULL,
                       populations = NULL,
                       data = NULL,
                       initialize = function(ldabase = LDABase$new(), granges = GRanges(), populations = NA){
                         self$ldabase <- ldabase
                         self$granges <- granges
                         self$populations <- populations
                       },set_data = function(){
                         self$data <- import1000GData(where = self$granges, which = self$populations, vcf_file = self$ldabase$vcf_file,
                                                      panel_file = self$ldabase$panel_file)
                         invisible(self$data)
                       }
                     )
)
