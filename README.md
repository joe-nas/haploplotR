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

In order to analyse linkage disequilibrium only in a window around our lead_snp we need to define a window around lead_snps. 
Here we take the lead_snps identifiers and create a GRanges object from it and subsequently resize the genomic position to a window.
For compatibility reason we set the seqlevelsStyle to "UCSC"
```r
lead_snps_gr<-GRanges(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37,lead_snps))
analysis_intervals_gr<-resize(lead_snps_gr,fix='center',width=6.5e5)
seqlevelsStyle(lead_snps_gr)<-"UCSC"
seqlevelsStyle(analysis_intervals_gr)<-"UCSC"
```

Next we put our ingredients into the LDABase, LDAImport and LDAnalysis1 objects.
Here for each of the analysis intervals and populations which are defined prior we import our data using the LDAImport object and by applying the $new() method. 
Subsequently, LDAnalysis1$new() is called. 
Here, we input the vector of lead_snps, define an rsquared cutoff which defaults to .8 and use our data as imported by LDAImport.
We generate results using the LDAnlaysis1$set_rsquared() and set_results() methods to calculate and filter which SNPs in the analysis windows are in linkage disequilibrium with the lead_snps. by meeting the cutoff value.
Here the set_rsquared() method calculates all rsquared values for lead snps with all other snps in the analysis window.
set_results() filters the rsquared data to meet the cutoff value.
```r

ldabase <- LDABase$new()
# workflow for only one interval/ single locus
dat <- LDAImport$new(ldabase = ldabase, granges = analysis_intervals_gr[1], 
                     populations = populations)$set_data()
res <- LDAnalysis1$new(lead_snps = lead_snps, cutoff = 0.8, lda_import = dat)
res$set_rsquared()
res$set_results()
```

The results can then be accessed by using res$results. res$results holds a Granges object containing the positions of  tag_snps annotated with:

- lead_snp_name
- tag_snp_name
- rsquared
- population
- lead_snp_pos



## To do:
- import1000GData.R needs a rewrite using VariantAnnotation as WhopGenome is no longer maintained
- write a decent howto, how to visualize
- tidy up roxygen
- tidy up remove: 
  - obsolete/analysis.R
  - obsolete/liloanalysis.R
  - obsolete/testscript.R
  - obsolete/LDAnalysis_class.R

