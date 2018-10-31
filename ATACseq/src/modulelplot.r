library(pheatmap)
pdf(file="module_plot.pdf",width=20)
load("../../ATACseq/results/peakdata_ATAC.Rdata")
load("../../ATACseq/results/peakdata_RNA.Rdata")
emat <- log2(as.matrix(bed)+1 )
rownames(emat) <- genename

load("../../ATACseq/results/modulelist.Rdata")

library(GSVA)

ord <- order(TISSUE,TYPE, AGE)
emat <- emat[,ord]

es <- gsva(emat, module)
## Estimating GSVA scores for 28 gene sets.
## Computing observed enrichment scores
## Estimating ECDFs with Gaussian kernels
## Using parallel with 8 cores
## 

anno <- data.frame(age=AGE[ord], strain=TYPE[ord],tissue=TISSUE[ord])
rownames(anno) <- colnames(es)
pheatmap(es, annotation_col=anno, cluster_cols=FALSE,show_colnames=F,main="RNAseq data")



dev.off()
