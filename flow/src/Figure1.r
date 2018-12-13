library(pheatmap)
library(reshape2)
library(ggplot2)
library(edgeR)
library(grid)
library(gridExtra)
library(gtable)

drawtable <- function(tab,title)
{
  t1 <- tableGrob(tab)
title <- textGrob(title)
padding <- unit(5,"mm")
  table <- gtable_add_rows(t1, heights = grobHeight(title) + padding,pos = 0)

table <- gtable_add_grob(table, title,  1, 1, 1, ncol(table))
  return(table)
}

source("../../ATACseq/src/library_function.r")
source("../../ATACseq/src/PVCA.r")


PCA <- function(emat,color.var,title)
  {
emat=emat[rowSums(emat)>0,]

pca=prcomp(t(emat),center=T,scale=T)

d <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)

xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])


dat=data.frame(PC1=as.numeric(pca$x[,1]),PC2=as.numeric(pca$x[,2]),tissue=color.var[,"TISSUE"],strain=color.var[,"STRAIN"],age=color.var[,"AGE"],gender=color.var[,"GENDER"])

i=3

p <- ggplot(dat, aes(PC1,PC2))+ geom_point(aes(color=dat[,i]))+labs(x=xl,y=yl,color=colnames(dat)[i],title=title)+ scale_colour_manual(  values=c("PBL"="red","spleen"= "blue", "memory"="green","naive"="yellow","BM"="orange"))
return(p)

  }





load("../../flow/data/flow_data_table.Rdata")

pdf(file="../results/Figure1.pdf",width=12)
flow[flow=="SPL"]="spleen"
flow[,4]=round(as.numeric(flow[,4])*7/30)
flow[flow[,4]==13,4]=12
flow[flow[,4]==17,4]=18
flow=flow[flow[,"tissue"]!="BM",]

### get rid of problematic middle age samples
flow=flow[flow[,"age"]!=6 ,]

meta=flow[,1:4]

tab=table(meta[,1],as.numeric(meta[,4]),meta[,3],meta[,2])

flowtable=grid.arrange(drawtable(tab[,,"F","B6"],"B6, F"),drawtable(tab[,,"M","B6"],"B6, M"),drawtable(tab[,,"F","NZO"],"NZO, F"),drawtable(tab[,,"M","NZO"],"NZO, M"),ncol=4)

### short table version ####
## tab=table(meta[,1],meta[,3],meta[,2])

## total <- function(x) {y=rbind(x,colSums(x));rownames(y)[nrow(y)]="total";return(y)}
## flowtable=grid.arrange(drawtable(total(tab[,,"B6"]),"B6"),drawtable(total(tab[,,"NZO"]),"NZO"),ncol=2)
####

colnames(meta)=c("TISSUE","STRAIN","GENDER","AGE")
dat=as.matrix(flow[,-(1:5)])
dat=matrix(as.numeric(dat),nr=nrow(dat),nc=ncol(dat))
colnames(dat)=colnames(flow)[-(1:5)]

dat=dat[,-c(26:36,48:51)]

res1=PVCA(t(dat[rowSums(is.na(dat))==0,]),meta[rowSums(is.na(dat))==0,],0.9)

p1=PCA(t(dat[rowSums(is.na(dat))==0,]),meta[rowSums(is.na(dat))==0,],"(C) PCA plot of each datset
flow")



#### ATAC data #####

load("../../ATACseq/data/ATACseqData.Rdata")

### normalization #####

y <- DGEList(counts=bed[,-(1:4)])
Y <- calcNormFactors(y,method="TMM")

STRAIN=TYPE

bed= cpm(Y, normalized.lib.sizes=TRUE)

source("../../ATACseq/src/remove_problemsample.r")

color.var=cbind(TISSUE,STRAIN,AGE,GENDER)

res2=PVCA(exp(bed[,]),color.var,0.8)

p2=PCA(bed,color.var,"
ATAC")

tab=table(TISSUE,as.numeric(AGE),GENDER,STRAIN)

ATACtable=grid.arrange(drawtable(tab[,,"F","B6"],"B6, F"),drawtable(tab[,,"M","B6"],"B6, M"),drawtable(tab[,,"F","NZO"],"NZO, F"),drawtable(tab[,,"M","NZO"],"NZO, M"),ncol=4)


load("../../RNAseq/data/RNAseqData.Rdata")

#### normalization #####

y <- DGEList(counts=bed)
Y <- calcNormFactors(y,method="TMM")


STRAIN=TYPE

bed= cpm(Y, normalized.lib.sizes=TRUE)

source("../../ATACseq/src/remove_problemsample.r")


color.var=cbind(TISSUE,STRAIN,AGE,GENDER)

res3=PVCA(exp(bed),color.var,1)

p3=PCA(bed,color.var,"
RNA")


tab=table(TISSUE,as.numeric(AGE),GENDER,STRAIN)

RNAtable=grid.arrange(drawtable(tab[,,"F","B6"],"B6, F"),drawtable(tab[,,"M","B6"],"B6, M"),drawtable(tab[,,"F","NZO"],"NZO, F"),drawtable(tab[,,"M","NZO"],"NZO, M"),ncol=4)


#####
factor=c(names(res1),names(res2),names(res3))
factor[factor=="resid"]=" Unexplained"
factor[factor=="GENDER"]="SEX"
ggdf=data.frame(data=rep(c("flow","ATAC","RNA"),each=length(res1)),value=c(res1,res2,res3),Factor=factor)
ggdf$data=factor(ggdf$data,levels=c("flow","ATAC","RNA"),ordered=TRUE)
#ggdf$Factor=factor(ggdf$Factor,levels=c("AGE","SEX","STRAIN","TISSUE"," Unexplained"))

p0=ggplot(ggdf, aes(x = data, y = value, fill = factor)) +theme_bw()+      geom_bar(stat = "identity") +theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Data")+ylab("proportion of variance")+ scale_fill_manual(name="",breaks = c("AGE", "SEX", "STRAIN","TISSUE"," Unexplained"),  values=c("grey","orange", "blue", "green","yellow"))+ggtitle("(B) proportion of variance explained by each factor")+ theme(plot.title = element_text(hjust= 0.5))
                                        #   +  scale_fill_manual(values = c("red","blue","green"))                                                                        

multiplot(NA,p1,NA,p2,p0,p3,cols=3)




dev.off()


pdf(file="table.pdf",width=10)
grid.arrange(textGrob("(A) number of flow samples                                                                             ",gp=gpar(fontsize=20)),flowtable,textGrob("(B) number of ATAC samples                                                                           ",gp=gpar(fontsize=20)),ATACtable,textGrob("(C) number of RNA samples                                                                           ",gp=gpar(fontsize=20)),RNAtable,ncol=1,heights=c(2,13,2,20,2,20))
dev.off()
