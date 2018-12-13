library(pheatmap)
library(edgeR)
library(ggplot2)

source("../../ATACseq/src/library_function.r")
source("../../ATACseq/src/PVCA.r")


total.number.plot=bmpathimmune=spleenpathimmune=pblpathimmune=memorypathimmune= naivepathimmune=bmpathwiki=spleenpathwiki=pblpathwiki= memorypathwiki= naivepathwiki=bmpathkegg=spleenpathkegg=pblpathkegg= memorypathkegg= naivepathkegg=vector("list",10)

balloonplot=vector("list",4)
YLIMMAX=6000

for(tid in 3:10)
  {
    source("../../ATACseq/src/matching_ATAC_RNA.r")

#    source("../../ATACseq/src/generate_pmat_fcmat.r")
    
    source("../../ATACseq/src/find_differential_peak_or_gene_strain.r")
  }
