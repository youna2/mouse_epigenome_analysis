load("../../RNAseq/data/RNAseqData.Rdata")
load("../../RNAseq/results/pmat_fcmat.Rdata")
EntrezRNA=annotation[,"Entrez.ID"]
annotationRNA=annotation
p.mat.RNA=p.mat
fc.mat.RNA=fc.mat
load("../../ATACseq/data/ATACseqData.Rdata")
load("../../ATACseq/results/pmat_fcmat.Rdata")
EntrezATAC=annotation[,"Entrez.ID"]
annotationATAC=annotation
p.mat.ATAC=p.mat
fc.mat.ATAC=fc.mat

common=intersect(unique(EntrezRNA),unique(EntrezATAC))
length(common)
matching.index=matrix(NA,nr=length(common),nc=2)
argmin <- function(x) return(which(x==min(x))[1])

for(i in 1:length(common))
  {
    matching.index[i,2]=match(common[i],EntrezRNA)
    
    temp=which(EntrezATAC %in% common[i])
    x=argmin(abs(as.numeric(annotationATAC[temp,"Distance.to.TSS"])))
    matching.index[i,1]=temp[x]
  }
mean(EntrezATAC[matching.index[,1]]==EntrezRNA[matching.index[,2]])
save(matching.index,file="../../ATACseq/results/matching_index.Rdata")

for(i in 1:ncol(fc.mat.ATAC))
  {

    plot(fc.mat.ATAC[matching.index[,1],i],fc.mat.RNA[matching.index[,2],i],main=colnames(fc.mat.RNA)[i],xlab="log FC ATAC",ylab="log FC RNA",pch='.')

    temp=cor.test(fc.mat.ATAC[matching.index[,1],i],fc.mat.RNA[matching.index[,2],i])
    legend("topright",legend=paste("cor=",signif(temp$estimate,2),"\n","p=",signif(temp$p.value,2)),box.lty=0)
           
  }
