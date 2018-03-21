require(R6)
LDAanalysis1 <- R6Class("LDAanalysis1",
                        public = list(
                          lda_import = NULL,
                          lead_snps = NULL,
                          cutoff = NULL,
                          rsquared = NULL,
                          results = NULL,
                          gr = NULL,
                          initialize = function(lead_snps = NA, cutoff = 0.8, lda_import = NA, gr = NULL){
                            self$lda_import <- lda_import
                            self$lead_snps <- lead_snps
                            self$cutoff <- cutoff
                            self$gr <- gr
                          },
                          set_rsquared = function(){
                            cat(paste0("Calculating R squared.\n"))
                            self$rsquared <- llply(self$lda_import, function(y){
                              calculateLD(lead_snps = lead_snps, xSnpMatrix = y$genotype)})
                            self$rsquared <<- Filter(function(x) {!is.null(x)}, self$rsquared)
                            invisible(self$rsquared)
                          },
                          set_results = function(){
                            cat(paste0("Finding HighLD intervals.\n"))
                            #self$rsquared <<- Filter(function(x) {!is.null(x)}, self$rsquared)
                            self$results <- Reduce(c, llply(names(self$rsquared), function(y){
                              identifyHighLD(rsquared = self$rsquared[[y]], info = self$lda_import[[y]]$info, cutoff = self$cutoff,
                                             population = y, gr = self$gr)}))
                            invisible(self$results)
                          }
                        ))
