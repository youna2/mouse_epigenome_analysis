heatmap <- function(x,title)
  {
    MAX=ceiling(quantile(abs(x),0.99))
    breaksList = seq(-MAX, MAX, by = 1)
    x[x>MAX]=MAX
    x[x< -MAX]= -MAX
    pheatmap(x,scale="none",annotation_col=anno,cluster_cols = FALSE,main=title,color = colorRampPalette(c("blue","white","red"))(length(breaksList)),breaks = breaksList,fontsize=12)
}


library(pheatmap)
load("../../ATACseq/results/FC_ATAC_RNA_NZO.Rdata")
ATAC.FC.NZO=ATAC.FC
RNA.FC.NZO=RNA.FC

load("../../ATACseq/results/FC_ATAC_RNA_B6.Rdata")
ATAC.FC.B6=ATAC.FC
RNA.FC.B6=RNA.FC


ATAC=cbind(ATAC.FC.B6,ATAC.FC.NZO,RNA.FC.B6,RNA.FC.NZO)
rownames(ATAC)=convert.mouse.entrez.to.human.symbol(mouseentrez)


anno <- data.frame(strain=as.matrix(data.frame(strsplit(colnames(ATAC)," ")))[3,],datatype=c(rep("ATAC",ncol(ATAC.FC.B6)+ncol(ATAC.FC.NZO)),rep("RNA",ncol(RNA.FC.B6)+ncol(RNA.FC.NZO))))

colnames(ATAC)=paste(colnames(ATAC),anno$datatype)
rownames(anno) <- colnames(ATAC)


load("../../ATACseq/results/inflammation_gene_list.Rdata")

inflammation3=as.matrix(read.csv("../../ATACseq/results/inflammation_gene_list.csv"))
myeloid.innate.immunity=as.matrix(read.csv("../../ATACseq/results/myeloid_innate_immunity_gene_list.csv"))

pdf(file="../../ATACseq/results/immune_gene_agign_change.pdf",width=20,height=34)
heatmap(ATAC[toupper(rownames(ATAC))%in% toupper(inflammation1),],"Inflammation I")

heatmap(ATAC[toupper(rownames(ATAC))%in% toupper(inflammation2),],"Inflammation II")

heatmap(ATAC[toupper(rownames(ATAC))%in% toupper(inflammation3),],"Inflammation genes from Kyungin")

heatmap(ATAC[toupper(rownames(ATAC))%in% toupper(myeloid.innate.immunity),],"Myeloid Innate Immunity genes from Kyungin")

dev.off()



