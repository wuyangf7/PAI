# PAI (promoter-anchored chromatin interactions)
This repository provides the scripts to implement the SMR and HEIDI methods to predict the promoter-anchored chromatin interactions using summary-level molecular quantitative trait locus (QTL) data. The package includes a shiny visualization tool to visualize and download the predicted interaction results computed using summary-level DNA methylation QTL data. 

## Reference

Wu Y, Qi T, Wang H, Zhang F, Zheng Z, Phillips-Cremins JE, Deary IJ, McRae AF, Wray NR, Zeng J, Yang J. Promoter-anchored chromatin interactions predicted by genetic analysis of epigenomic data. BioRxiv.

https://www.biorxiv.org/content/10.1101/580993v1

## Getting Started

Our analytical approach relies on two recently developed methods, i.e., the summary-data–based Mendelian randomization (SMR) test and the test for heterogeneity in dependent instruments (HEIDI). To install SMR, you can download the smr_Linux.zip package from the SMR website (http://cnsgenomics.com/software/smr/#Overview), which contains a standalone (i.e., statically linked) 64-bit Linux executable file smr_Linux. We strongly recommend using this static executable because it is well-optimized and no further installation is required. 

Here we provide an example to predict the promoter-anchored chromatin interactions using mQTL summary data from peripheral blood samples. In this case, we will need the mQTL summary data in BESD format (http://cnsgenomics.com/software/smr/#DataManagement). We then focus on testing for pleiotropic associations of a DNAm site in the promoter region of a gene with all the other DNAm sites within 2 Mb of the focal promoter in either direction. 

```
smr --bfile mydata --beqtl-summary myDNAm --extract-exposure-probe DNAmInPromoter --beqtl-summary myDNAm  --out myPAIs
```

* --bfile  reads individual-level SNP genotype data (in PLINK binary format) from a reference sample for LD estimation, i.e. .bed, .bim, and .fam files.
* --beqtl-summary the first one reads mQTL summary data as the exposure and the second one reads mQTL summary data from as the outcome. 
* --extract-exposure-probe extracts a subset of exposure DNAm probes in promoter regions for PAI analysis.
