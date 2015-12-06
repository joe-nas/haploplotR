

ycoordFUN <- function(txs){
  unlist(lapply(1:length(txs), function(i){
    sum(unlist(lapply(seq.int(i),function(j){
      overlap(txs[i],txs[j])
    })))
  }))
}

overlap <- function(r1, r2){
  start(r1) <= end(r2) && start(r2) <= end(r1)
}



