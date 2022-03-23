# haploplotR:  Project to visualize linkage disequilibrium from 1000 genomes data

![rs2235371](https://user-images.githubusercontent.com/4353093/159236477-4170a802-f476-4a18-bb08-a036a98b7ada.png)

Example plot using the IRF6 locus. Red bars in the population panels represent tag snps, 
snps wich are in linkage disequilibrium (rsquared >= .8) with lead snps indicated by blue bars.
The p63BS panel indicates p63 binding sites which overlap with lead or tag snps (red) in at least one population 
and therefore are at least in part in linkage disequilibrium with the lead snp. Or as indicated by black bars are not in linkage disequilibrium with a lead or tag snp.

For this particular locus lead snps are snps indicated by gwas focusing on non syndromic cleft lip with or without cleft palate (NSCL/P). 
Given that gwas only indicates a genomic region, further study of that region needs to be done. By analysing linkage disequilibrium in that region
and bringing it into context with p63, a tf regulating facial development one might find regulatory regions like enhancers being involved in NSCL/P.

## requirements:
- 1000 genomes vcf files: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- 1000 genomes panel: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
- lead_snps: vector of dbsnp identifiers e.g. c("rs79997038","rs79084855","rs560426","rs66515264")
- population: vector of population identifiers as in panel file e.g. c("CEU", "HAN", "ACB")

## How to:

First a few dependencies and useful objects need to be loaded
```r {cmd}
#loading dependencies
library(devtools)
library(GenomicRanges)
library(grid)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)
devtools::load_all("path_to/haploplotR")
```

Next we specify 1000genomes vcf files and panel file paths. Furthermore target snps and populations are defined.
```r
# specify 1000genomes vcf files and panel file paths
vcffile<-"path_to/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
panelfile<-"path_to/integrated_call_samples_v3.20130502.ALL.panel"
 
# create vectors defining lead_snps as well as 1000genomes population identifiers

lead_snps<-c("rs79997038","rs79084855","rs560426")
populations <- c("CEU",'CHB',"LWK")
```
in order to analyse linkage disequilibrium only in a window around our lead_snp we need to define a window around lead_snps.
Here we take the lead_snps identifiers and create a GRanges object from it and subsequently resize the genomic position to a window.
For compatibility reason we set the seqlevelsStyle to "UCSC"
```r

lead_snps_gr<-GRanges(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37,lead_snps))
analysis_intervals_gr<-resize(lead_snps_gr,fix='center',width=6.5e5)
seqlevelsStyle(lead_snps_gr)<-"UCSC"
seqlevelsStyle(analysis_intervals_gr)<-"UCSC"

```


Generally speaking you define a GRanges object describing the region of interest, 
1000genomes populations you are interested in eg. c("CEU", "HAN", "ACB") and 
finally a vector containing lead-snps c("rs79997038","rs79084855","rs560426","rs66515264").
Lead-snps are the snps you start calculation pairwise linkage disequilibrium from. 

## To do:
- import1000GData.R needs a rewrite using VariantAnnotation as WhopGenome is no longer maintained
- write a decent howto
- tidy up roxygen
- tidy up remove: 
  - obsolete/analysis.R
  - obsolete/liloanalysis.R
  - obsolete/testscript.R
  - obsolete/LDAnalysis_class.R

