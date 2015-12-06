collapseTranscriptsBy <- function(txs, by){
  txs <- GenomicRanges::as.data.frame(txs)
  txs$symbol <- select(Homo.sapiens::Homo.sapiens, keys = as.character(txs$tx_name),
                       keytype = "TXNAME", columns = "SYMBOL")[,"SYMBOL"]


  txs_s <- Filter(nrow, split(txs, list(txs$"type", do.call("$",list(txs, by)))))
  lapply(txs_s, function(x){
    tmp <- with(x, reduce(GRanges(Rle(seqnames),
                                  IRanges(start, end),
                                  strand = unique(x$strand))))
    GRanges(Rle(seqnames(tmp)), ranges(tmp), strand = unique(x$strand),
            type = rep(unique(x$type), length(tmp)),
            symbol = rep(unique(x$symbol), length(tmp))
    )
  })
}
