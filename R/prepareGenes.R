prepareGenes <- function(where){
  df <- GenomicRanges::as.data.frame(GenomicRanges::sort(biovizBase::crunch(TxDb.Hsapiens.UCSC.hg19.knownGene,where)))
  df$SYMBOL <- select(Homo.sapiens, keys = as.character(df$gene_id), keytype = "GENEID",
                      columns = "SYMBOL")[,"SYMBOL"]
  stepping <- setYCoordsFUN(Reduce(rbind,lapply(split(df,df$SYMBOL), function(x){
    data.frame(seqnames = unique(x$seqnames),start = min(x$start)-20000,end = max(x$end)+20000,
               SYMBOL = unique(x$SYMBOL), strand = unique(x$strand)[1])})))
  df <- merge(df,stepping[c("SYMBOL","y")], by="SYMBOL",all.x = T)
  df <- df[complete.cases(df),]
  genes_df <- GenomicRanges::as.data.frame(Reduce(c,lapply(split(df, list(df$type, df$SYMBOL)),
                  function(x){
                    tmp <- reduce(with(x, IRanges(start = start,end = end)))
                    with(x, GRanges(Rle(unique(seqnames)),
                                    IRanges(start(tmp), end(tmp)),
                                    strand = unique(x$strand),
                                    SYMBOL = rep(unique(x$SYMBOL), length(tmp)),
                                    type = rep(unique(x$type), length(tmp)),
                                    y = rep(unique(x$y), length(tmp))))
                  })))
  list(genes = genes_df, genes_ranges = with(stepping,
                                             GRanges(Rle(seqnames),
                                                     IRanges(start = start+20000, end = end-20000),
                                                     strand = strand, SYMBOL = SYMBOL,
                                                     y = y)))
}
