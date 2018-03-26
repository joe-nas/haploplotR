identifyHighLD <- function(rsquared, info, gr = NULL, cutoff = 0.8, population = NULL){
  # index in rsquared
  highld_idx <- which(rsquared >= cutoff, arr.ind = T)
  tag_snp_idx <- highld_idx[,1] # corresponds with info as well
  lead_snp_idx <- highld_idx[,2] # corresponds wtih rsquared only

  # index in info
  # ifelse tests used to circumvent problems ocrruring when generating interfals generated on idx 1 the end of the range
  tag_interval_start_v <- ifelse(tag_snp_idx > 1, tag_snp_idx - 1, tag_snp_idx) # start position of tag_interval
  tag_interval_end_v <- ifelse(tag_snp_idx < nrow(rsquared), tag_snp_idx + 1, tag_snp_idx) # stop position of tag_interval


  # compiling vectors for results data.frame
  #meta columns
  lead_snp_names_v <- colnames(rsquared)[highld_idx[,2]]
  #lead_snp_names_v <- names(which(rsquared >= 0.8, arr.ind=T)[,2])
  tag_snp_names_v <- info[highld_idx[,1],]$ids
  rs_value_v <- rsquared[highld_idx]
  # lead_snp_pos_v <- info$pos[which(lead_snp_names_v == info$ids)] ## THIS NEEDS TO BE DEBUGGED
  tag_snp_pos_v <- info$pos[highld_idx[,1]]

  #DEBUB
  # print(tag_snp_names_v)
  # print(lead_snp_names_v)

  meta_df <- data.frame(lead_snp = lead_snp_names_v, tag_snp = tag_snp_names_v, rsquared = rs_value_v,
                        population = population, tag_snp_pos = tag_snp_pos_v)

  if(!is.null(gr) && length(gr)>0){
    a_interval_str <- with(gr, paste(seqnames,start,end,sep=":"))
    meta_df <- data.frame(meta_df, a_interval = a_interval_str)
  }
  cat("This is identifyHighLD.R \n")

  #GRanges
  chromosomes <- paste("chr", info$chromosome_num[tag_snp_idx], sep = "")
  start <- info$pos[tag_interval_start_v]
  end <- info$pos[tag_interval_end_v]

  GRanges(Rle(chromosomes), IRanges(start,end), mcols = meta_df)
}
