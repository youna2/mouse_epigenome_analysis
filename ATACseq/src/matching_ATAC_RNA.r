### This code is old, and the object bed is not included in the ATACseqData2.Rdata or RNAseqData2.Rdata anymore. If to use, must fix this


library(edgeR)
setwd("../results")
previous.dir=getwd()
utissue=c( "spleen", "BM", "CD8.memory","CD8.naive", "PBL")         

p.cutoff=0.05

if(tid>=5)
  {
    commonpattern=TRUE
    if(tid==5) int.or.slope="B6 vs NZO common"
    if(tid==6) int.or.slope="B6 vs NZO opposite"
    if(tid==7) int.or.slope="only NZO significant"
    if(tid==8) int.or.slope="only B6 significant"
    if(tid==9) int.or.slope="within B6"
    if(tid==10) int.or.slope="within NZO"
    if(tid==11) int.or.slope="all"
    topgene=int.or.slope

    if(tid==5 | tid>=9)
      {
        positivepeak="opening"
        negativepeak="closing"
      }else{
        if(tid==8)
          {
            positivepeak="higher in B6"
            negativepeak="lower in B6"
          }else{
            positivepeak="higher in NZO"
            negativepeak="lower in NZO"
          }
      }
    
  }else{
    commonpattern=FALSE
    if(tid==3)
      {
        typeinteraction=TRUE
        int.or.slope="slope"
        topgene="B6_vs_NZO_slope"
      
      }else{
        typeinteraction=FALSE
        int.or.slope="intercept"
        topgene="B6_vs_NZO_intercept"
      }
    positivepeak="higher in NZO"
    negativepeak="lower in NZO"

  }


###############################

load("../../RNAseq/results/RNAseqData2.Rdata")

load(paste("../../RNAseq/results/pmat_fcmat_",int.or.slope,".Rdata",sep=""))

p.mat.RNA=p.mat.strain
fc.mat.RNA=fc.mat.strain

EntrezRNA=annotation[,"Entrez.ID"]
annotationRNA=annotation


load("../../ATACseq/results/ATACseqData2.Rdata")

## y <- DGEList(counts=bed[,-(1:4)])
## Y <- calcNormFactors(y,method="TMM")
## bed <- cpm(Y, normalized.lib.sizes=TRUE)

load(paste("../../ATACseq/results/pmat_fcmat_",int.or.slope,".Rdata",sep="")) 

p.mat.ATAC=p.mat.strain
fc.mat.ATAC=fc.mat.strain

EntrezATAC=annotation[,"Entrez.ID"]
annotationATAC=annotation
bedATAC=bed

common=intersect(unique(EntrezRNA),unique(EntrezATAC))
length(common)
matching.index=matrix(NA,nr=length(common),nc=2)
argmin <- function(x) return(which(x==min(x))[1])

for(j in 1:length(common))
  {
    matching.index[j,2]=match(common[j],EntrezRNA)
    
    temp=which(EntrezATAC %in% common[j])
    x=argmin(abs(as.numeric(annotationATAC[temp,"Distance.to.TSS"])))
    matching.index[j,1]=temp[x]
  }
mean(EntrezATAC[matching.index[,1]]==EntrezRNA[matching.index[,2]])

dim(fc.mat.ATAC[matching.index[,1],])
dim(p.mat.ATAC[matching.index[,1],])
dim(annotationATAC[matching.index[,1],])

dim(fc.mat.RNA[matching.index[,2],])
dim(p.mat.RNA[matching.index[,2],])
dim(annotationRNA[matching.index[,2],])

commoncolumn=intersect(colnames(p.mat.RNA),colnames(p.mat.ATAC))
p.mat.RNA=p.mat.RNA[,match(commoncolumn,colnames(p.mat.RNA))]
fc.mat.RNA=fc.mat.RNA[,match(commoncolumn,colnames(fc.mat.RNA))]
p.mat.ATAC=p.mat.ATAC[,match(commoncolumn,colnames(p.mat.ATAC))]
fc.mat.ATAC=fc.mat.ATAC[,match(commoncolumn,colnames(fc.mat.ATAC))]


p.mat.strain=p.mat.ATAC[matching.index[,1],]
fc.mat.strain=fc.mat.ATAC[matching.index[,1],]
annotation=annotationATAC[matching.index[,1],]
bed=bedATAC[matching.index[,1],]

sum(p.mat.strain<p.cutoff)
p.mat.strain[fc.mat.strain * fc.mat.RNA[matching.index[,2],]<0]=1
#p.mat.strain[p.mat.RNA[matching.index[,2],]>p.cutoff]=1
sum(p.mat.strain<p.cutoff)

###
convert.to.human.gene.name=F
if(convert.to.human.gene.name)
  {

row.bed=convert.mouse.entrez.to.human.symbol(annotation[ ,"Entrez.ID"])

bed=bed[!is.na(row.bed),]
genename=row.bed[!is.na(row.bed)]


save(bed,AGE,TYPE,TISSUE,genename,file="peakdata.Rdata")
}
