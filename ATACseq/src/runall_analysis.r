library(pheatmap)
library(edgeR)
library(ggplot2)
library(preprocessCore)

for(BTID in 2:5)
  {
    
load(paste("../../ATACseq/data/ATACseqData",BTID,".Rdata",sep=""))
source("../../ATACseq/src/library_function.r")
source("../../ATACseq/src/PVCA.r")
### normalization #####


librarysize=colSums(as.matrix(bed[,-(1:4)]))
if(BTID<6)
  {
    bed[,-(1:4)] <- normalize.quantiles(as.matrix(bed[,-(1:4)]))
  }else{
    tmp=unique(TISSUE)
    for(i in 1:length(tmp))
      {
        bed[,-(1:4)][,TISSUE==tmp[i]] <- normalize.quantiles(as.matrix(bed[,-(1:4)][,TISSUE==tmp[i]]))
      }
  }
y <- DGEList(counts=bed[,-(1:4)])
Y <- calcNormFactors(y,method="none") ##method="TMM"
#librarysize=y$samples[,"lib.size"]


pdf(file=paste("../results",BTID,"/PCA_ATAC.pdf",sep=""))
STRAIN=TYPE
color.var=cbind(TISSUE,STRAIN,AGE,GENDER)

if(BTID==6)
{
  res=PVCA(bed[,-(1:4)],color.var,0.8)
  PlotPVCA(res, "")
}
bed= cpm(Y, normalized.lib.sizes=TRUE)
PCA(bed,color.var)

bed.noTMM= cpm(Y, normalized.lib.sizes=TRUE)
bed.noTMM=bed.noTMM[,order(librarysize)]
librarysize=librarysize[order(librarysize)]
## par(mfrow=c(4,1))
## for(i in seq(from=1,to=ncol(bed.noTMM),length.out=10)) hist(bed.noTMM[,i][bed.noTMM[,i]<=quantile(bed.noTMM[,i],0.99)],breaks=1000,main=paste("lib=",signif(librarysize[i],2),"var=",round(var(bed.noTMM[,i]),2)))


top.bot=top.bottom(bed.noTMM,0.9,0.1)
par(mfrow=c(2,2))
plot(AGE,librarysize,xlab="age",ylab="library size")#age vs library size

plot(top.bot[,1],librarysize,xlab="average read counts of top 90%",ylab="library size")
temp=cor.test(top.bot[,1],librarysize)
legend("topleft",legend=c(paste("r=",signif(temp$estimate,2)),paste("p=",signif(temp$p.value,2))),box.lty=0)

plot(top.bot[,2],librarysize,xlab="average read counts of bottom 10%",ylab="library size")
temp=cor.test(top.bot[,2],librarysize)
legend("topleft",legend=c(paste("r=",signif(temp$estimate,2)),paste("p=",signif(temp$p.value,2))),box.lty=0)

#plot(colSums(bed.noTMM),librarysize)

plot(apply(bed.noTMM,2,var),librarysize,xlab="varaince",ylab="library size")
temp=cor.test(apply(bed.noTMM,2,var),librarysize)
legend("topleft",legend=c(paste("r=",signif(temp$estimate,2)),paste("p=",signif(temp$p.value,2))),box.lty=0)

#plot(apply(bed.noTMM,2,median),librarysize)

dev.off()

LIBRARYSIZE=log(y$samples[,"lib.size"])
DistanceToTSS= 2500  ### 10^7
if(BTID<7)
  {
    source("../../ATACseq/src/remove_problemsample.r")
  } else{
y0forheatmap=bed
bed=log(bed+min(bed[bed>0]))
selectgene=rep(TRUE,nrow(bed))
annotation.orig=annotation
save(bed,LIBRARYSIZE,annotation.orig,annotation,AGE,GENDER,TISSUE,STRAIN,TYPE,SAMPLEMOUSEID,file=paste("../results",BTID,"/ATACseqData2.Rdata",sep=""))
tid=9
source("../../ATACseq/src/generate_pmat_fcmat.r")

  }

## Do timeseries analysis separately for each tissue type and gender and find peaks that are increasing/decreasing with age and do pathway enrichment analysis for the identified genes. There are four options :
## tid = 1 : do time series analysis in B6 strain. In each tissue type and gender, fit y~ c0+ age*c1  where y is chromatin accessibility of a peak and find peaks with significant p-values for the coefficient c1. 
## tid = 2 : do time series analysis in NZO strain. In each tissue type and gender, fit y~ c2 + age*c3  where y is chromatin accessibility of a peak and find peaks with significant p-values for the coefficient c3.
## tid = 3 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c0 is different from c2 and identify peaks where c0 and c2 are significantly different.
## tid = 4 : In each tissue type, for each peak, fit y~ c0+ age*c1 for B6 and fit y~ c2 + age*c3 for NZO and test if c1 is different from c3 and identify peaks where c1 and c3 are significantly different.
## tid = 5 : In each tissue type, for each peak, fit y~ c0+ c_strain+ c_gender+ age*c1 and find peaks with significant p-values for the coefficient c1. This is to find common aging pattern across gender and strain. 

YLIMMAX=10000
if(BTID<7) source("../../ATACseq/src/runall_analysis_second.r")
}
source("../../ATACseq/src/merge_results.r")
