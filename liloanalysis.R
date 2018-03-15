library(devtools)
library(GenomicRanges)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(Homo.sapiens)
library(doMC)
library(plyr)

if(Sys.getenv("LOGNAME") == "jfalck"){
  devtools::load_all("/scratch/jfalck/git/haploplotR")
  vcffile<-"/scratch/jfalck/1000G/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  panelfile<-"/scratch/jfalck/1000G/integrated_call_samples_v3.20130502.ALL.panel"
  registerDoMC(cores=16)
  lead_snps<-c("rs79997038","rs79084855","rs560426","rs66515264","rs751399","rs10498466","rs12569773","rs4615961","rs10983654",
               "rs1512262","rs11119348","rs12566152","rs2865509","rs11698025","rs59138205","rs59626211","rs227727","rs227731",
               "rs9904526","rs9891446","rs112924906","rs643310","rs28474857","rs10886036","rs1898349","rs7017665","rs13266917",
               "rs139971355","rs143478693","rs188674070","rs12895971","rs1047057","rs1649202","rs1877432","rs2808672","rs7543674",
               "rs11119438","rs11699548","rs2144129","rs62291848","rs28809250","rs16957824","rs7208324","rs6503169","rs2177787",
               "rs116790447","rs77826974","rs61076166","rs10985356","rs61873036","rs72728755","rs17242358","rs10818050","rs7362194",
               "rs7679350","rs59308021","rs72809921","rs72809924","rs12028175","rs4920339","rs357537","rs987525","rs642961",
               "rs17085106","rs7590268","rs10863790","rs7078160","rs9574565","rs1258763","rs17760296","rs2294426","rs861020",
               "rs13041247","rs742071","rs7632427","rs12543318","rs8001641","rs1873147","rs8076457","rs4441471","rs1373453",
               "rs7846606","rs7820074","rs4132699","rs13542","rs813218","rs3815854","rs4703516","rs7950069","rs5765956",
               "rs1536895")[c(8:11,28)]

  populations<-c("ACB","ASW","ESN","GWD","LWK","MSL","YRI","CLM","MXL","PEL","PUR","CDX","CHB","CHS","JPT","KHV","CEU",
                 "FIN","GBR","IBS","TSI","BEB", "GIH","ITU","PJL","STU")
}else{
  devtools::load_all("/home/SSD-Data/Projects/haploplotR")
  vcffile<-"/home/SSD-Data/1000Genomes/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
  panelfile<-"integrated_call_samples_v3.20130502.ALL.panel"
  registerDoMC(cores=4)
  lead_snps<-c("rs79997038","rs79084855","rs560426","rs66515264","rs751399","rs10498466","rs12569773","rs4615961","rs10983654",
               "rs1512262","rs11119348","rs12566152","rs2865509","rs11698025","rs59138205","rs59626211","rs227727","rs227731",
               "rs9904526","rs9891446","rs112924906","rs643310","rs28474857","rs10886036","rs1898349","rs7017665","rs13266917",
               "rs139971355","rs143478693","rs188674070","rs12895971","rs1047057","rs1649202","rs1877432","rs2808672","rs7543674",
               "rs11119438","rs11699548","rs2144129","rs62291848","rs28809250","rs16957824","rs7208324","rs6503169","rs2177787",
               "rs116790447","rs77826974","rs61076166","rs10985356","rs61873036","rs72728755","rs17242358","rs10818050","rs7362194",
               "rs7679350","rs59308021","rs72809921","rs72809924","rs12028175","rs4920339","rs357537","rs987525","rs642961",
               "rs17085106","rs7590268","rs10863790","rs7078160","rs9574565","rs1258763","rs17760296","rs2294426","rs861020",
               "rs13041247","rs742071","rs7632427","rs12543318","rs8001641","rs1873147","rs8076457","rs4441471","rs1373453",
               "rs7846606","rs7820074","rs4132699","rs13542","rs813218","rs3815854","rs4703516","rs7950069","rs5765956",
               "rs1536895")[c(8:11,28)]

  populations<-c("ACB","ASW","ESN","GWD","LWK","MSL","YRI","CLM","MXL","PEL","PUR","CDX","CHB","CHS","JPT","KHV","CEU",
                 "FIN","GBR","IBS","TSI","BEB", "GIH","ITU","PJL","STU")
  populations <- populations[populations %in% c("CEU",'CHB',"LWK")]

}





lead_snps_gr<-GRanges(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37,lead_snps))
analysis_intervals_gr<-resize(lead_snps_gr,fix='center',width=6.5e5)
seqlevelsStyle(lead_snps_gr)<-"UCSC"
seqlevelsStyle(analysis_intervals_gr)<-"UCSC"
genome(lead_snps_gr)<-NA
genome(analysis_intervals_gr)<-NA

TGData<-llply(analysis_intervals_gr,function(x){
  gr<-x#locationwherersuarediscalculated
  data<-import1000GData(where=x,which=populations,vcf_file=vcffile,panel_file=panelfile)#1000Gdata
  rsquared<-Filter(function(y){!is.null(y)},llply(data,function(dat) calculateLD(lead_snps = lead_snps, dat$genotype)))#rsforgivenpopulationandlead_snps
  # highld <- Reduce(c,llply(names(rsquared),function(y){
  #  identifyHighLD(rsquared=rsquared[[y]],data[[y]]$info,population=y)
  # }))
  list(gr=gr,data=data,rsquared=rsquared, lead_snps=lead_snps)#,
       # highld = highld)#returnallthedataaslist
},.parallel=TRUE)

saveRDS(TGData, file="TGData.rds")

results<-llply(TGData,function(x){
  Reduce(c,llply(names(x$rsquared),function(y){
    identifyHighLD(rsquared=x$rsquared[[y]],x$data[[y]]$info,population=y)
  }))
})


results<-with(TGData[[1]],Reduce(c,llply(names(rsquared),function(y){
    identifyHighLD(rsquared=rsquared[[y]],data[[y]]$info,population=y)})))
