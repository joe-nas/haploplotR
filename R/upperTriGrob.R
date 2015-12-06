#' grid Grob representing the upper triangle of pairwise SNP statistics like R^2
#'
#' @param r2 Symmetric numeric matrix with range [0,1] representing pairwise R^2 values.
#' @param r2d Depth of shown results. r2d = 100 specifies for each SNP only the
#' adjacent 100 pairwise statistics should be shown.
#' @param mapping DataFrame with columns 'pos' and 'ids' specifying SNP positions
#' and SNP names.
#' @param lead Vector of SNP identifiers for which tag SNPs are identifierd and
#' which are be highlighted in the mapping panel.
#' @param r2_cutoff The cutoff value used to identify tag SNPs.
#' @param label A string which will be plotted in the upper left, e.g. 'CEU'.
#' @return gList
#' @examples
#' where <- GRanges(Rle('chr1'), IRanges(start = 209789729, end = 210189729)
#' mygenes <- genesGrob(where = where, what = TxDb.Hsapiens.UCSC.hg19.knownGene)
#' pushViewport(viewport)
#' grid.draw(mygenes)



upperTriGrob <- function(r2, r2d, mapping, r2_cutoff = 0.8, lead = NULL,
                         y = unit(5/6, 'npc'), xscale, label, colorscheme, ...){

  mapping <- prepareMapping(r2 = r2, lead = lead, r2_cutoff = r2_cutoff,
                            mapping = mapping)

  grid::gList(
    grid::rasterGrob(prepareRaster(mat = r2, depths = r2d, colorscheme = colorscheme),
                     vp = haploplot.vp(y = y), interpolate = F),

    grid::textGrob(
      label = label , x = 0, y = 1, just  = 'left',
      gp = grid::gpar(fontface = 'bold'),
      vp = grid::viewport(y = y, height = unit(1/6, 'npc'),
                          width = unit(sqrt(2), 'snpc')), ...),

    grid::segmentsGrob(
      x0 = unit(sort(mapping$pos), 'native'),
      x1 = unit(sort(1/length(mapping$pos)) * seq.int(length(mapping$pos)),'npc'),
      y0 = unit(0.6, 'native'), y1 = unit(1, 'native'),
      gp = gpar(alpha = mapping$alpha, col = mapping$color),
      vp = mapping.vp(x = mapping$pos, y = y, just = 'top', xscale = xscale, ...)),

    grid::segmentsGrob(
      x0 = unit(sort(mapping$pos), 'native'),
      x1 = unit(sort(mapping$pos), 'native'),
      y0 = unit(0.2, 'native'), y1 = unit(0.6, 'native'),
      gp = gpar(alpha = mapping$alpha, col = mapping$color),
      vp = mapping.vp(x = mapping$pos, y = y, just = 'top', xscale = xscale, ...))

    )
}
