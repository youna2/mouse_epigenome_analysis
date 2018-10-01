library(pheatmap)
library(edgeR)
library(ggplot2)

load("../../ATACseq/data/ATACseqData.Rdata")
source("../../ATACseq/src/library_function.r")
source("../../ATACseq/src/PVCA.r")
### normalization #####

y <- DGEList(counts=bed[,-(1:4)])
Y <- calcNormFactors(y,method="TMM")


pdf(file="../results/PCA_ATAC.pdf")
STRAIN=TYPE
color.var=cbind(TISSUE,STRAIN,AGE,GENDER)

res=PVCA(bed[,-(1:4)],color.var,0.8)
PlotPVCA(res, "")

bed[,-(1:4)]= cpm(Y, normalized.lib.sizes=TRUE)

PCA(bed[,-(1:4)],color.var)
dev.off()


## Do timeseries analysis separately for each tissue type and gender and find peaks that are increasing/decreasing with age and do pathway enrichment analysis for the identified genes. There are four options :
## tid = 1 : do time series analysis in B6 strain. In each tissue type and gender, fit y~ c0+ age*c1  where y is chromatin accessibility of a peak and find peaks with significant p-values for the coefficient c1. 
## tid = 2 : do time series analysis in NZO strain. In each tissue type and gender, fit y~ c2 + age*c3  where y is chromatin accessibility of a peak and find peaks with significant p-values for the coefficient c3.
## tid = 3 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c0 is different from c2 and identify peaks where c0 and c2 are significantly different.
## tid = 4 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c1 is different from c3 and identify peaks where c1 and c3 are significantly different.
## tid = 5 : In each tissue type, for each peak, fit y~ c0+ c_strain+ c_gender+ age*c1 and find peaks with significant p-values for the coefficient c1. This is to find common aging pattern across gender and strain. 



tid=1

if(tid==1) selB6=TRUE
if(tid==2) selB6=FALSE

YLIM=c(-30000,30000)
if(tid>2) source("../../ATACseq/src/find_differential_peak_or_gene_strain.r") else source("../../ATACseq/src/find_differential_peak_or_gene.r")


