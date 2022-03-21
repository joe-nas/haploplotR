# haploplotR:  Project to visualize linkage disequilibrium from 1000 genomes data

## requirements:
- 1000 genomes vcf files
- 1000 genomes panel file
- lead_snps: vector of dbsnp identifiers e.g. c("rs79997038","rs79084855","rs560426","rs66515264")
- population: vector of population identifiers as in panel file e.g. c("CEU", "HAN", "ACB")

Generally speaking you define a GRanges object describing the region of interest, 
1000genomes populations you are interested in eg. c("CEU", "HAN", "ACB") and 
finally a vector containing lead-snps c("rs79997038","rs79084855","rs560426","rs66515264").
Lead-snps are the snps you start calculation pairwise linkage disequilibrium from. 
