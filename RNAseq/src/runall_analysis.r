library(pheatmap)
library(ggplot2)
library(edgeR)

load("../../RNAseq/data/RNAseqData.Rdata")
source("../../ATACseq/src/library_function.r")
source("../../ATACseq/src/PVCA.r")
#### normalization #####

y <- DGEList(counts=bed)
Y <- calcNormFactors(y,method="TMM")


pdf(file="../results/PCA_RNA.pdf")
STRAIN=TYPE
color.var=cbind(TISSUE,STRAIN,AGE,GENDER)

res=PVCA(bed,color.var,1)
PlotPVCA(res, "")

bed= cpm(Y, normalized.lib.sizes=TRUE)
PCA(bed,color.var)
dev.off()


## Do timeseries analysis separately for each tissue type and gender and find genes that are increasing/decreasing with age and do pathway enrichment analysis for the identified genes. There are four options :
## tid = 1 : do time series analysis in B6 strain. In each tissue type and gender, fit y~ c0+ age*c1  where y is chromatin accessibility of a peak and find genes with significant p-values for the coefficient c1. 
## tid = 2 : do time series analysis in NZO strain. In each tissue type and gender, fit y~ c2 + age*c3  where y is chromatin accessibility of a peak and find genes with significant p-values for the coefficient c3.
## tid = 3 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c0 is different from c2 and identify genes where c0 and c2 are significantly different.
## tid = 4 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c1 is different from c3 and identify genes where c1 and c3 are significantly different.

YLIMMAX=4000
source("../../ATACseq/src/runall_analysis_second.r")


