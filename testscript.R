library(devtools)
library(GenomicRanges)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
devtools::load_all('/home/jonas/R/haploplotR/')


# library(org.Hs.eg.db)
# interesting.anno <- addGeneIDs(annotatedPeak=interesting.anno,
#                                orgAnn="org.Hs.eg.db",
#                                IDs2Add="symbol")
#



target_gr <- GRanges(Rle('chr1'), IRanges(start = 209789729, end = 210189729))
vcffile <- "ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
panelfile <- "integrated_call_samples_v3.20130502.ALL.panel"


target_data <- get_SnpMatrix_from_vcf(where = target_gr, which = "CEU",
                                      vcf_file = vcffile,
                                      panel_file = panelfile)

rsquared_ceu <- snpStats::ld(target_data$genotype, target_data$genotype, stats = 'R.squared')


colscheme <- colorRampPalette(c('#fff7fb','#3690c0','#023858'))(101)



populations <- c('CEU','CHB','JPT','YRI')
names(populations) <- c('CEU','CHB','JPT','YRI')
target_dat <- lapply(populations, function(pop){
  get_SnpMatrix_from_vcf(where = target_gr, which = pop, vcf_file = vcffile,
                         panel_file = panelfile)
})


target_rsquared <- lapply(target_dat, function(dat){
  snpStats::ld(dat$genotype, dat$genotype, stats = 'R.squared')
})

my_mapping <- Reduce(rbind,lapply(seq_along(target_dat), function(i){
  prepareMapping(r2 = target_rsquared[[i]],mapping =  target_dat[[i]]$info,r2_cutoff = 0.8,
                 lead = c('rs10863790','rs10863790','rs11119348', 'rs12566152','rs861020','rs10863790'))
}))

my_mapping <- my_mapping[!duplicated(my_mapping),]
my_mapping <- my_mapping[order(my_mapping$pos),]


grid::grid.newpage()
lapply(seq_along(populations), function(i){
  grid::grid.draw(
    upperTriGrob(r2 = target_rsquared[[populations[i]]], r2d = 250, mapping = target_dat[[populations[i]]]$info,
                 xscale = c(start(target_gr), end(target_gr)), r2_cutoff = 0.8,
                 label = populations[i], colorscheme = colscheme,y = unit((6-i)/6, 'npc'),
                 lead = c('rs10863790','rs10863790','rs11119348', 'rs12566152','rs861020','rs10863790'))
  )
})


grid.draw(gTree(vp = viewport(height = unit(1/6, 'npc'),width = unit(sqrt(2), 'snpc'),
                              y = unit(1/6, 'npc'), xscale = c(start(target_gr), end(target_gr))),
                children = gList(xaxisGrob(at = round(seq(start(target_gr), end(target_gr), length.out = 12))[-c(1,12)],
                                           main = T, gp = gpar(fontsize = 10)),
                                 segmentsGrob(x0 = unit(my_mapping$pos, 'native'), x1 = unit(my_mapping$pos, 'native'),
                                              y0 = unit(1.15, 'native'), y1 = unit(1.2, 'native'),
                                              gp  = gpar(col=my_mapping$color, fill = my_mapping$color,
                                                         alpha = my_mapping$alpha)),
                                 genesGrob(where = target_gr, what =  TxDb.Hsapiens.UCSC.hg19.knownGene))))






gf <- frameGrob()
gf <- packGrob(gf, gTree(vp = viewport(height = unit(1/6, 'npc'),width = unit(sqrt(2), 'snpc'),
                                       y = unit(1/6, 'npc'), xscale = c(start(target_gr), end(target_gr))),
                         children = gList(xaxisGrob(at = round(seq(start(target_gr), end(target_gr), length.out = 12))[-c(1,12)],
                                                    main = T, gp = gpar(fontsize = 10)),
                                          segmentsGrob(x0 = unit(my_mapping$pos, 'native'), x1 = unit(my_mapping$pos, 'native'),
                                                       y0 = unit(1.15, 'native'), y1 = unit(1.2, 'native'),
                                                       gp  = gpar(col=my_mapping$color, fill = my_mapping$color,
                                                                  alpha = my_mapping$alpha)),
                                          genesGrob(where = target_gr, what =  TxDb.Hsapiens.UCSC.hg19.knownGene))))

grid.draw(gf)


target_data$genotype
head(target_data$info)










## --------- Analysis ---------
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(plyr)
library(doMC)
registerDoMC(2)


lead_snps <- read.table("~/../Storage/Documents/p63_latex/section_results/tables/lead_snps.csv",sep = ',', header = T, stringsAsFactors = F)
lead_snps <- snpsById(SNPlocs.Hsapiens.dbSNP142.GRCh37, lead_snps$lead.SNPs)


lead_snps_res <- resize(lead_snps, fix = 'center', width = 4e5)
lead_snps_red <- reduce(lead_snps_res)
lead_snps_red <- resize(lead_snps_red, fix = 'center', width = 4e5)
seqlevelsStyle(lead_snps_red) <- "UCSC"
seqlevelsStyle(lead_snps) <- "UCSC"
genome(lead_snps_red) <- NA
genome(lead_snps) <- NA


setwd("~/../Storage/Downloads/1000G/")

haploplot <- function(where, leads){
  where <- GRanges(Rle(as.character(seqnames(where))), IRanges(start(where), end(where)))
  colscheme <- colorRampPalette(c('#fff7fb','#3690c0','#023858'))(101)
  vcffile <- "ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  panelfile <- "integrated_call_samples_v3.20130502.ALL.panel"

  populations <- c('CEU','CHB','JPT','YRI')
  names(populations) <- c('CEU','CHB','JPT','YRI')

  target_dat <- lapply(populations, function(pop){
    get_SnpMatrix_from_vcf(where = where, which = pop, vcf_file = vcffile,
                           panel_file = panelfile)
  })

  target_rsquared <- lapply(target_dat, function(dat){
    snpStats::ld(dat$genotype, dat$genotype, stats = 'R.squared')
  })

  hld <- Reduce(c, lapply(populations, function(pop){
    hldfun(rs = target_rsquared[[pop]], info = target_dat[[pop]]$info,
           snp_names = leads, pop = pop, cutoff = 0.8,
           chr = paste("chr", target_dat[[pop]]$info$chromosome_num[1], sep = ""))
  }))

  my_mapping <- Reduce(rbind,lapply(seq_along(target_dat), function(i){
    prepareMapping(r2 = target_rsquared[[i]],
                   mapping =  target_dat[[i]]$info,r2_cutoff = 0.8,
                   lead = leads)
  }))

  my_mapping <- my_mapping[!duplicated(my_mapping),]
  my_mapping <- my_mapping[order(my_mapping$pos),]
  my_mapping <- my_mapping[my_mapping$col != "black",]

  grid::grid.newpage()
  grid.draw(textGrob(label = with(where, paste(seqnames,paste(start, end, sep = "-"),sep = ":")),
                     x = 0.5, y = 1, just = "top", gp = gpar(fontsize = 18)))

  lapply(seq_along(populations), function(i){
    grid::grid.draw(
      upperTriGrob(r2 = target_rsquared[[populations[i]]], r2d = 200,
                   mapping = target_dat[[populations[i]]]$info,
                   xscale = c(start(where), end(where)), r2_cutoff = 0.8,
                   label = populations[i], colorscheme = colscheme,
                   y = unit((6-i)/6, 'npc'),
                   lead = leads)
    )
  })


  gf <- frameGrob()
  genesgrob <-genesGrob(where = where, what =  TxDb.Hsapiens.UCSC.hg19.knownGene)
  summarygrob <- segmentsGrob(x0 = unit(my_mapping$pos, 'native'),
                              x1 = unit(my_mapping$pos, 'native'),
                              y0 = unit(1.175, 'native'),
                              y1 = unit(1.225, 'native'),
                              gp  = gpar(col = my_mapping$color,
                                         fill = my_mapping$color,
                                         alpha = my_mapping$alpha))
  chipldgrob <- chipLDGrob(chip = p63_hg19, LD = hld, where = where)

  gf <- packGrob(gf, gTree(vp = viewport(height = unit(1/6, 'npc'),
                                         width = unit(sqrt(2), 'snpc'),
                                         y = unit(1/6, 'npc'),
                                         xscale = c(start(where), end(where))),
                           children = gList(xaxisGrob(at = round(seq(start(where), end(where), length.out = 12))[-c(1,12)],
                                                      main = T, gp = gpar(fontsize = 10)),
                                            genesgrob, summarygrob, chipldgrob)
                           ))
  grid.draw(gf)


}


genesGrob(where = lead_snps_red, what =  TxDb.Hsapiens.UCSC.hg19.knownGene)

library(Cairo)
# l_ply(1:33,function(i){
#   Cairo(file = paste(paste("region",i , sep = "_") ,"png", sep = '.'),
l_ply(seq_along(lead_snps$RefSNP_id),function(i){
  Cairo(file = paste(lead_snps$RefSNP_id[i],"png", sep = '.'),
        width = 789, height = 493, units = "px", type = 'png', bg = 'white')
  grid.newpage()
  capture.output(try(haploplot(resize(lead_snps[i], width = 400000, fix = "center"),
                               lead_snps$RefSNP_id)), file = NULL)
  dev.off()
  gc(verbose = F)
},.parallel = F, .progress = "text")


grid.newpage()
try(haploplot(lead_snps_red[5],lead_snps$RefSNP_id))


hldfun <- function(rs, info, snp_names, chr, pop, cutoff = 0.8){
  rows <- rownames(rs) %in% snp_names
  diffs <- diff(info[,1])
  centers <- info[,1]
  starts <- c(centers-c(0,diffs))
  ends <- c(centers+c(diffs,0))
  if(sum(rows) == 1){
    subspace <- rs[rows,]
    idx <- which(subspace >= cutoff, arr.ind = T, useNames = F)
    lead_name <- info[rows,"ids"]
    tag_name <- names(subspace)[idx]
    rs_val <- subspace[idx]
    meta <- data.frame(tag = tag_name, lead = lead_name,
                       rs = rs_val, population = pop,
                       stringsAsFactors = F)
    gr <- GRanges(Rle(chr), IRanges(starts[idx], ends[idx]))
    mcols(gr) <- meta
    gr
  }else if(sum(rows) >= 2){
    subspace <- rs[rows,]
    idx <- which(subspace >= cutoff, arr.ind = T, useNames = F)
    lead_name <- rownames(subspace[idx[,1],])
    tag_name <- colnames(subspace[,idx[,2]])
    rs_val <- subspace[idx]
    meta <- data.frame(tag = tag_name, lead = lead_name,
                       rs = rs_val, population = pop,
                       stringsAsFactors = F)

    pos <- info[info$ids %in% tag_name, "pos"]
    gr <- GRanges(Rle(chr), IRanges(starts[idx[,2]],ends[idx[,2]]))
    mcols(gr) <- meta
    gr
  }else{
      GRanges()
    }
}

chipLDGrob <- function(chip, LD, where){
  chip$highlight <- "black"
  ov <- findOverlaps(query = chip, subject =  LD)
  chip$highlight[unique(queryHits(ov))] <- "red"
  hl <- as.data.frame(subsetByOverlaps(chip,where))
  rects <- with(hl, rectGrob(x = unit(start,"native"),
                    width = unit(end-start,"native"),
                    y = unit(1.1, 'native'),
                    height = unit(0.05, 'native'),
                    just = "bottom",
                    gp = gpar(col = highlight,
                              fill = highlight)))
  annos <- textGrob(label = rownames(hl),
                    x = unit(with(hl,(start+end)/2), 'native'),
                    y = unit(1.075, 'native'),
                    just = 'top',
                    gp = gpar(fontsize = 8))
  gList(rects, annos)
}


test <- debugfun(where = lead_snps_red[5],leads = lead_snps$RefSNP_id)
highldregresions <- Reduce(c,lapply(c("CEU","CHB","JPT","YRI"), function(pop){
  hldfun(rs = test$tr$CEU, info = test$td[[pop]]$info, snp_names = lead_snps$RefSNP_id,
         chr = paste("chr",test$td[[pop]]$info$chromosome_num[1], sep = ""),pop = pop)
}))

chipLDGrob(chip = p63_hg19, LD = highldregresions, where = lead_snps_red[1])
grid.newpage(grid.draw(rectGrob(x = unit(start(bla),"native") , width =  unit(width(bla),"native"),
                   y = 0.5,height = 1, gp = gpar(col = bla$highlight, fill = bla$highlight) , vp = viewport(xscale = c(start(lead_snps_red[1]),
                                            end(lead_snps_red[1]))))))



library(TxDb.Hsapiens.UCSC.hg19.knownGene)
prepareGenes(with(lead_snps_red[23], GRanges(seqnames,IRanges(start,end))))
grid.newpage()
grid.draw(
  genesGrob(what = TxDb.Hsapiens.UCSC.hg19.knownGene,
            where = with(lead_snps_red[23], GRanges(seqnames,IRanges(start,end))))
)


prepareGenes(target_gr)

## --------------- lead snps and regions

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
lead_snps
seqlevelsStyle(lead_snps) <- 'UCSC'
genome(lead_snps) <- NULL
lead_snps <- annotatePeakInBatch(lead_snps, AnnotationData = TSS.human.GRCh37)

library(org.Hs.eg.db)
lead_snps.anno <- addGeneIDs(annotatedPeak=lead_snps,
                             orgAnn="org.Hs.eg.db",
                             IDs2Add="symbol")

lead_snps.appendix <- as.data.frame(lead_snps.anno)
head(lead_snps.appendix)
rownames(lead_snps.appendix) <- NULL
lead_snps.appendix
write.table(as.data.frame(lead_snps.appendix), "/home/Storage/Documents/p63_latex/section_appendix/tables/appendix.lead_snps_hg19.tab", quote = F, sep = '\t',row.names = F, col.names = T)



## -------------------------------------

lead_snps <- resize(lead_snps, fix = 'center', width = 4e5)
names(lead_snps) <- lead_snps$RefSNP_id
seqlevelsStyle(lead_snps) <- 'UCSC'



populations <- c('CEU','CHB','JPT','YRI')
names(populations) <- c('CEU','CHB','JPT','YRI')


vcffile <- "ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
panelfile <- "integrated_call_samples_v3.20130502.ALL.panel"



genotypes <- llply(lead_snps, function(lead){
  llply(populations, function(pop){
    get_SnpMatrix_from_vcf(where = lead, which = pop, vcf_file = vcffile,
                           panel_file = panelfile)
  }, .parallel = T)
})





highLD <- function(rs, gt, lead, population, cutoff = 0.8){
  #lead <- unique(leads[lead %in% rownames(rs)])
  lead_gr <- lead
  lead <- lead$RefSNP_id
  lead <- unique(lead[lead %in% rownames(rs)])
  #tags <- unique(which(rs[leads,] >= cutoff, arr.ind = T)[,2]
  if(!length(lead) == 0){
    print("in if")
    tags <- unique(which(rs[lead,] >= cutoff, arr.ind = T))
    tags <- if(length(lead) > 1){
      tags[,2]
    }else{
      print("in else")
      tags
    }

    diffs <- diff(gt$info[,1])
    centers <- gt$info[,1]
    starts <- c(centers-c(0,diffs))
    ends <- c(centers+c(diffs,0))

    gr <- GRanges(Rle(as.character(seqnames(lead_gr))), IRanges(starts, ends),
                  lead_snp = names(lead_gr),tag_snp = gt$info[,2], ld = rs[lead,],
                  population = population)
    names(gr) <- gt$info[,2]
    gr[tags]
  }else{
    print("granges empty")
    GRanges()
  }
}





highld <- llply(lead_snps, function(snp){
  gc(verbose = F)
  llply(names(genotypes[[names(snp)]]),function(p){
    gc(verbose = F)
    highLD(rs = with(genotypes[[names(snp)]][[p]],
                     snpStats::ld(genotype, genotype, stats = 'R.squared')),
           gt = genotypes[[names(snp)]][[p]],
           lead = snp, population = p, cutoff = 0.8)
  }, .parallel = T)
}, .progress = 'text')


highld <- Reduce(c,unlist(highld))







library(rtracklayer)
source("../p63/scripts.R")
p63_1 <- with(read.table("../p63/GSM439776_p63_1_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4, p63BS = paste("p63s1",1:length(V2),sep="_")))
p63_2 <- with(read.table("../p63/GSM439777_p63_2_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4, p63BS = paste("p63s2",1:length(V2),sep="_")))
p63_3 <- with(read.table("../p63/GSM538930_p63_3_peaks.bed.gz",header = F),GRanges(Rle(V1), IRanges(start=V2,end=V3),score = V4, p63BS = paste("p63s3",1:length(V2),sep="_")))
p63 <- reduce(c(p63_1,p63_2,p63_3))

ch <- import.chain(ChainFile("../p63/hg18ToHg19.over.chain"))
p63_hg19 <- unlist(liftOver(x = p63,chain = ch))
p63_hg19$p63BS <- 1:length(p63_hg19)
names(p63_hg19) <- 1:length(p63_hg19)

write.table(as.data.frame(p63_hg19), "/home/Storage/Documents/p63_latex/section_appendix/tables/appendix.p63BS_hg19.tab", quote = F, sep = '\t',row.names = F, col.names = T)





motifs <- subseqfun(p63_hg19)
motifs$name <- paste('motif',1:length(motifs),sep="")
motifs$motif <- "motif"

write.table(as.data.frame(motifs)[-7], "/home/Storage/Documents/p63_latex/section_appendix/tables/appendix.p63motifs_hg19.tab", quote = F, sep = '\t',row.names = F, col.names = T)


motifs_anno <- annoGR(motifs)

library(EnsDb.Hsapiens.v75)
library(ChIPpeakAnno)
EnsDb.Hsapiens.v75
myanno <- annoGR(EnsDb.Hsapiens.v75, feature = 'gene')



## -------------
ol <- findOverlaps(p63_hg19, highld)
interesting <- GRanges(p63_hg19[queryHits(ol)],
                       meta = cbind(elementMetadata(p63_hg19[queryHits(ol)]),
                                    elementMetadata(highld[subjectHits(ol)])))


H3k4me1 <- with(read.table("../histone/wgEncodeBroadHistoneNhekH3k4me1StdPk.txt",header = F),
                GRanges(Rle(V2),IRanges(V3,V4), histone = rep("H3k4me1",length(V1))))
H3k4me3 <- with(read.table("../histone/wgEncodeBroadHistoneNhekH3k4me3StdPk.txt",header = F),
                GRanges(Rle(V2),IRanges(V3,V4), histone = rep("H3k4me3",length(V1))))
H3k27ac <- with(read.table("../histone/wgEncodeBroadHistoneNhekH3k27acStdPk.txt",header = F),
                GRanges(Rle(V2),IRanges(V3,V4), histone = rep("H3k27ac",length(V1))))



customAnno <- function(gr1, gr2, label){
  gr1$label = F
  ol <- findOverlaps(gr1, gr2)
  GRanges(gr1[queryHits(ol)], meta = cbind(elementMetadata(gr1[queryHits(ol)]),
                                           elementMetadata(gr2[subjectHits(ol)])))
}




interesting$motif <- "no"
ol <- findOverlaps(interesting, motifs)
interesting[queryHits(ol)]$motif <- 'yes'


interesting$H3k4me1 <- "no"
ol <- findOverlaps(interesting, H3k4me1)
interesting[queryHits(ol)]$H3k4me1 <- 'yes'


interesting$H3k4me3 <- "no"
ol <- findOverlaps(interesting, H3k4me3)
interesting[queryHits(ol)]$H3k4me3 <- 'yes'


interesting$H3k27ac <- "no"
ol <- findOverlaps(interesting, H3k27ac)
interesting[queryHits(ol)]$H3k27ac <- 'yes'


names(interesting) <- NULL




data(TSS.human.GRCh37)
library(EnsDb.Hsapiens.v75)
anno <- annoGR(EnsDb.Hsapiens.v75, 'gene')
interesting.anno <- annotatePeakInBatch(interesting, AnnotationData=anno)

library(org.Hs.eg.db)
interesting.anno <- addGeneIDs(annotatedPeak=interesting.anno,
                               orgAnn="org.Hs.eg.db",
                               IDs2Add="symbol")

## -------------- appendix
res_tab <- as.data.frame(interesting.anno)
res_tab



## -------------- text + format
res_text <- res_tab[,c(1,2,3,7,10,11,12,13,14,16,24)]

rownames(res_text) <- NULL
res_text$p63BS <- with(res_text, paste(seqnames,":",start,"-",end,sep = ""))
res_text <- res_text[c(1,2,12,11,10,4,5,6,7,8,9)]
res_text <- res_text[!duplicated(res_text),]
res_text$seqnames <- as.integer(gsub("chr","",res_text$seqnames))
res_text <- res_text[order(res_text$seqnames, res_text$start,res_text$meta.lead_snp,res_text$meta.population),]

res_text<-res_text %>% spread(meta.population,meta.population)
res_text$population <- with(res_text, gsub(". NA","",gsub("NA, ","",paste(CEU,CHB,JPT,YRI, sep=", "))))
res_text <- res_text[!colnames(res_text)%in%c('seqnames', 'start','CEU', 'CHB','JPT','YRI')]
res_text <- res_text[c(1,3,2,4,9,5,6,7,8)]
colnames(res_text) <- c("p63 BS","Ensembl ID", "Gene Symbol", "lead SNP", "Population","p63 Motif", "H3k4me1", "H3k4me3", "H3k27ac")

res_text$`p63 Motif` <- ifelse(res_text$`p63 Motif` == "yes", "+","-")
res_text$H3k4me1 <- ifelse(res_text$H3k4me1 == "yes", "+","-")
res_text$H3k4me3 <- ifelse(res_text$H3k4me3 == "yes", "+","-")
res_text$H3k27ac <- ifelse(res_text$H3k27ac == "yes", "+","-")
res_text

stargazer(res_text,summary = F, rownames = F, float.env = 'sidewaystable',font.size = "footnotesize")

library(stargazer)



## appendix + format

appendix <- res_tab

head(appendix)
rownames(appendix) <- NULL
appendix <- appendix[c(1,2,3,4,6,7,8,9,10,11,12,13,14,16,20,21,24)]
colnames(appendix) <- c("chr","start","end","width", "p63 BS","lead SNP","tag SNP","r squared","population", "motif", "H3k4me1","H3k4me3", "H3k27ac", "feature","insideFeature","distancetoFeature",'symbol')
write.table(appendix, "/home/Storage/Documents/p63_latex/section_appendix/tables/appendix.results.bs_lead_tag.tab", quote = F, sep = '\t',row.names = F, col.names = T)
head(appendix)



annotatePeakInBatch(interesting, AnnotationData = motifs_anno)

interesting_anno <- annotatePeakInBatch(interesting, AnnotationData = myanno)





library(Homo.sapiens)
columns(Homo.sapiens)
mytxs_ids <- select(Homo.sapiens, keys = as.character(mytxs$tx_name), keytype = "TXNAME", columns = "ENSEMBL")

library(biovizBase)
mytxs <- biovizBase::crunch(TxDb.Hsapiens.UCSC.hg19.knownGene, target_gr)
collapseTranscriptsBy <- function(txs, by){
  txs <- as.data.frame(txs)
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


res <- collapseTranscriptsBy(mytxs, "gene_id")



Reduce(unlist,res)
