#' Import 1000 Genomes vcf chunks and convert them to SnpMatrix objects
#'
#' @param where GRanges object specifying chromosome, start and end position of the target region.
#' @param which Target population(s). Something like "CEU" or c("CEU", "CHB").
#' @param vcf_file local file or url to vcf file.
#' @param panel_file local file or url to panel file.
#' @return Object of class SnpMatrix.
#' @examples
#' COMPLETE EXAMPLES


import1000GData <- function(where, which, vcf_file, panel_file, ...){
  # capture.output({
  panel_file <- read.table(panel_file, header = T, stringsAsFactors = F)
  panel_file_pops <- split(panel_file, panel_file$pop)[which]

  chromosome_str <- as.character(seqnames(where))
  chromosome_num <- gsub('chr', '', chromosome_str)
  start_pos <- start(where)
  end_pos <- end(where)

  gr_str <- paste(chromosome_str, start_pos, end_pos, sep = ":")
  sink("/dev/null")
  res <- llply(panel_file_pops, function(x){
    vcf_handle <- with(where, WhopGenome::vcf_open(sprintf(vcf_file, chromosome_str)))
    WhopGenome::vcf_setregion(vcffh = vcf_handle, chromosome_num, start_pos, end_pos)

    which <- x$sample

    #print(x)

    WhopGenome::vcf_selectsamples(vcf_handle, which)
    WhopGenome::vcf_getregion(vcf_handle)

    #print(WhopGenome::vcf_getregion(vcf_handle))
    ## create info df
    ids <- c()
    pos <- c()
    repeat{
      WhopGenome::vcf_parseNextSNP(vcf_handle)
      ids <- c(ids, WhopGenome::vcf_getID(vcf_handle))
      pos <- c(pos, WhopGenome::vcf_getPos(vcf_handle))
      if(is.null(WhopGenome::vcf_getPos(vcf_handle))){
        break
      }
    }
    info <- data.frame(pos, ids, chromosome_num, stringsAsFactors = F)

    ## get heterocygosity
    WhopGenome::vcf_restartregion(vcf_handle)
    ishet <- matrix(data=as.integer(0), nrow = length(which), ncol = length(pos),
                    dimnames = list(which,1:length(pos)))
    WhopGenome::VCF_snpmat_diplo_bial_ishet_unfiltered(vcf_handle, ishet)
    ishet <- ishet[,!colnames(ishet) %in% '-1']

    ## get alternative allele presence
    WhopGenome::vcf_restartregion(vcf_handle)
    hasalt <- matrix(data=as.integer(0), nrow = length(which), ncol = length(pos),
                     dimnames = list(which,1:length(pos)))
    WhopGenome::VCF_snpmat_diplo_bial_hasalt_unfiltered(vcf_handle, hasalt)
    hasalt <- hasalt[,!colnames(hasalt) %in% '-1']

    ## close file handler
    WhopGenome::vcf_close(vcf_handle)

    ## keep this columns
    keep <- intersect(colnames(hasalt), colnames(ishet))

    ## setting row/colnames
    info <- info[info[,1] %in%  keep,]
    rownames(info) <- NULL

    genotype <- haploplotR::toRawSnpMat(ishet[,keep], hasalt[,keep])
    dimnames(genotype) <- list(which, info$ids)

    ## create SnpMatrix object snpStats
    require(snpStats)
    genotype <- new('SnpMatrix', genotype)

    # })
    cat("This is import1000GData \n")
    list(genotype = genotype, info = info, gr = where, gr_str = gr_str)
  }, .parallel = T, warn.conflicts = FALSE)
  sink("/dev/stdout")
}
