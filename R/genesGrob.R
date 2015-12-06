#' grid Grob representing gene models
#'
#' @param where GRanges object specifying chromosome, start and end position of the target region.
#' @param what TxDb object, e.g.TxDb.Hsapiens.UCSC.hg19.knownGene
#' @return gList
#' @examples
#' where <- GRanges(Rle('chr1'), IRanges(start = 209789729, end = 210189729))
#' mygenes <- genesGrob(where = where, what = TxDb.Hsapiens.UCSC.hg19.knownGene)
#' pushViewport(viewport)
#' grid.draw(mygenes)



genesGrob <- function(where, what, ...){

  # functions as native coordinate system in viewports
  xrange <- with(where, c(start(where), end(where)))

  pg <- prepareGenes(where)
  mygenes <- pg[['genes']]
  mygenes_ranges <- pg[['genes_ranges']]

  myarrows <- as.data.frame(Reduce(f = rbind, x = lapply(
    seq_along(pg[['genes_ranges']]), function(i) {
      arrowFUN(where = where, narr = 80, gr = pg[['genes_ranges']][i])
    })
  ))


  exon <- with(mygenes[mygenes$type == 'exon',],
               grid::rectGrob(x = grid::unit(start, 'native'), y = grid::unit(y, 'native'),
                        height = grid::unit(0.5, 'char'), width = grid::unit(width, 'native'),
                        gp = grid::gpar(fill = 'black', col = 'black', lwd = 0.5),
                        vp = genes2.vp(xscale = xrange, yscale = c(0,max(mygenes$y)))
                        )
               )



  cds <- with(mygenes[mygenes$type == 'cds',],
              rectGrob(x = unit(start, 'native'),
                       y = unit(y, 'native'),
                       height = unit(0.5, 'char'),
                       width = unit(width, 'native'),
                       gp = gpar(fill = '#008837', col = '#008837', lwd = 0.5),
                       vp = genes2.vp(xscale = xrange,
                                      yscale = c(0,max(mygenes$y)))
                       )
              )


  utr <- with(mygenes[mygenes$type == 'utr',],
              rectGrob(x = unit(start, 'native'),
                       y = unit(y, 'native'),
                       height = unit(0.5, 'char'),
                       width = unit(width, 'native'),
                       gp = gpar(fill = '#7b3294', col = '#7b3294', lwd = 0.5),
                       vp = genes2.vp(xscale = xrange,
                                      yscale = c(0,max(mygenes$y)))
                       )
              )

  gaps <- with(mygenes[mygenes$type == 'gap',],
               segmentsGrob(x0 = unit(start, 'native'),
                            x1 = unit(end, 'native'),
                            y0 = unit(y, 'native'),
                            y1 = unit(y, 'native'),
                            gp = gpar(fill = 'black', col = 'black',lwd = 1),
                            vp = genes2.vp(xscale = xrange,
                                           yscale = c(0,max(mygenes$y)))
                            )
               )


  myarrows <- segmentsGrob(x0 = unit(myarrows$starts, 'native'),
                           x1 = unit(myarrows$ends, 'native'),
                           y0 = unit(myarrows$y, 'native'),
                           y1 = unit(myarrows$y, 'native'),
                           gp = gpar(fill = 'black', col = 'black', lwd = 1),
                           arrow = arrow(angle = 45,
                                         ends = 'first',
                                         type = 'open',
                                         length = unit(0.04, 'npc')),
                           vp = genes2.vp(xscale = xrange,
                                          yscale = c(0,max(mygenes$y)))
                           )


  anno <- textGrob(label = mygenes_ranges$SYMBOL,
                   check.overlap = F,
                   x = unit(start(mygenes_ranges)+ width(mygenes_ranges)/2, 'native'),
                   y = unit(mygenes_ranges$y + 0.5 , 'native'),
                   gp = gpar(fontsize = 8),
                   vp = genes2.vp(xscale = xrange,
                                  yscale = c(0,max(mygenes$y))))



  gList(gaps,
        myarrows,
        exon, cds, utr, anno)
}
