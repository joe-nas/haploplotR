#' internal functions used in genesGrob to produce arrows
#'
#' @param where GRanges object specifying chromosome, start and end position of the target region.
#' @param narr Number of arrows on xrange
#' @param gr GRanges object specifying the location of one transcript (start to
#' end) and the y position as produced by ycoordFUN
#' @return matrix A matrix comprising start, end, and y positions of arrows
#' @examples
#' where <- GRanges(Rle('chr1'), IRanges(start = 209789729, end = 210189729))



arrowFUN <- function(where, narr = 80, gr){

  step <- length(start(where):end(where))/narr
  arr <- seq(start(gr),end(gr), by=step)
  if(as.logical(strand(gr) == '+')){
    ends <- arr[-c(1,length(arr))]-2/step
    starts <- ends + 2
    y <- rep(gr$y, length(starts))
    return(cbind(starts, ends, y))
  }else{
    ends <- arr[-c(1,length(arr))]+2/step
    starts <- ends - 2
    y <- rep(gr$y, length(starts))
    return(cbind(starts, ends, y))
  }
}

