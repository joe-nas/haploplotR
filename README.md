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
Generally speaking you define a GRanges object describing the region of interest, 
1000genomes populations you are interested in eg. c("CEU", "HAN", "ACB") and 
finally a vector containing lead-snps c("rs79997038","rs79084855","rs560426","rs66515264").
Lead-snps are the snps you start calculation pairwise linkage disequilibrium from. 

## To do:
- write a decent howto
- tidy up roxygen
- tidy up remove: 
  - obsolete/analysis.R
  - obsolete/liloanalysis.R
  - obsolete/testscript.R
  - obsolete/LDAnalysis_class.R

