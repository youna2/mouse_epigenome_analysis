library(pheatmap)
library(reshape2)
library(ggplot2)


source("library_function.r")
source("../../ATACseq/src/PVCA.r")
load("../data/flow_data_table.Rdata")
flow[flow=="SPL"]="spleen"
flow[,4]=round(as.numeric(flow[,4])*7/30)
flow[flow[,4]==13,4]=12
flow[flow[,4]==17,4]=18
flow=flow[flow[,"tissue"]!="BM",]

### get rid of problematic middle age samples
flow=flow[flow[,"age"]!=6 ,]

meta=flow[,1:5]
colnames(meta)=c("TISSUE","STRAIN","GENDER","AGE","SAMPLEID")
dat=as.matrix(flow[,-(1:5)])
dat=matrix(as.numeric(dat),nr=nrow(dat),nc=ncol(dat))
colnames(dat)=colnames(flow)[-(1:5)]
res=PVCA(t(dat[rowSums(is.na(dat))==0,]),meta[rowSums(is.na(dat))==0,1:4],0.9)
pdf(file="../results/PVCA_flow.pdf")
PlotPVCA(res, "")
PCA(t(dat[rowSums(is.na(dat))==0,]),meta[rowSums(is.na(dat))==0,1:4])
dev.off()
#### 3. Start analysis ###################
colnames(dat)[26:36]
colnames(dat)[48:51]
dat=dat[,-c(26:36,48:51)]
colSums(is.na(dat))
bed=t(dat[rowSums(is.na(dat))==0,])
bed=bed[-(1:3),]
meta=meta[rowSums(is.na(dat))==0,]

TISSUE=meta[,1]
STRAIN=meta[,2]
GENDER=meta[,3]
AGE=as.numeric(meta[,4])
SAMPLEMOUSEID=meta[,5]
################
TRANSFORM=F

dat=as.data.frame(flow)

for(i in c(4,6:ncol(dat))) dat[,i]=as.numeric(as.matrix(dat[,i]))


dat$CD4.CD8.ratio=as.numeric(flow[,"CD4.Viable"])/as.numeric(flow[,"CD8.Viable"])
dat$CD4.CD8.Lymphocytes.ratio=as.numeric(flow[,"CD4.Lymphocytes"])/as.numeric(flow[,"CD8.Lymphocytes"])
#dat$CD4.Naive.EMRA=as.numeric(flow[,"CD4.Naive"])+as.numeric(flow[,"CD4.EMRA"])
#dat$CD8.Naive.EMRA=as.numeric(flow[,"CD8.Naive"])+as.numeric(flow[,"CD8.EMRA"])

dat=dat[!is.na(dat[,2]),]


### draw global plot of cell proportions ###
convert.numeric <- function(x)
  {
    x[x==3]="03"
    x[x==6]="06"
    return(x)
  }

pdf(file="../results/global_flow_plot_problem_samples_removed.pdf",width=15,height=15)

dat0=dat
dat0$new=paste(convert.numeric(dat0$age),as.character(dat0$Sample),sep="-")
dat0$new=factor(dat0$age)
dat0$B=dat0$B.Lymphocytes
dat0$CD4=dat0$CD4.Lymphocytes
dat0$CD8=dat0$CD8.Lymphocytes

dat0=dat0[dat0[,"tissue"]!="BM",]

## par(mfrow=c(2,2))
## plot(dat0$B.Viable,dat0$B.Lymphocytes,xlab="B.Viable",ylab="B.Lymphocytes")
## abline(0,1)

## plot(dat0$CD4.Viable,dat0$CD4.Lymphocytes,xlab="CD4.Viable",ylab="CD4.Lymphocytes")
## abline(0,1)
## plot(dat0$CD8.Viable,dat0$CD8.Lymphocytes,xlab="CD8.Viable",ylab="CD8.Lymphocytes")
## abline(0,1)

## plot(dat0$CD4.CD8.ratio,dat0$CD4.CD8.Lymphocytes.ratio,xlab="CD4.Viable/CD8.Viable",ylab="CD4.Lymphocytes/CD8.Lymphocytes")
## abline(0,1)


#ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("B.Viable","CD4.Viable","CD8.Viable"))
#subplot(ggdf)


ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("B","CD4","CD8"))
p1=subplot(ggdf,"(A) percentage of sorted T cells")


#ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("CD4.CD8.ratio","CD4.CD8.Lymphocytes.ratio"))
#subplot3(ggdf,"NA")



ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "CD4.Naive", "CD4.CM", "CD4.EM", "CD4.EMRA"  ))
p2=subplot(ggdf,"")

#ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c(  "CD4.CM", "CD4.EM" ))
#subplot(ggdf)
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "CD8.Naive", "CD8.CM", "CD8.EM", "CD8.EMRA"  ))
p3=subplot(ggdf,"")    



#ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "X..Viable","Single.Cells.Viable.Cells.Granulocytes.Fre.Viable.Cells" , "Single.Cells.Viable.Cells.No.Granulocytes.Monocytes.Fre.Parent"   , "Single.Cells.Viable.Cells.No.Granulocytes.Monocytes.Fre.Viable.Cells"   ))
#subplot3(ggdf, "Single.Cells.Viable.Cells.")


#ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c(  "CD8.CM", "CD8.EM"  ))
#subplot(ggdf)

ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("IL7R.B","IL7R.CD4","IL7R.CD8"))
p4=subplot4(ggdf,"IL7R","(B) percentage of IL7R+ cells in sorted T cells" )
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "IL7R.CD4.Naive", "IL7R.CD4.CM", "IL7R.CD4.EM", "IL7R.CD4.EMRA"  ))
p5=subplot4(ggdf,"IL7R","")
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "IL7R.CD8.Naive", "IL7R.CD8.CM", "IL7R.CD8.EM", "IL7R.CD8.EMRA"  ))
p6=subplot4(ggdf,"IL7R","")
    
## ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("CD38.B","CD38.CD4","CD38.CD8"))
## subplot2(ggdf,"CD38")
    
## ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "CD38.CD4.Naive", "CD38.CD4.CM", "CD38.CD4.EM", "CD38.CD4.EMRA"  ))
## subplot2(ggdf,"CD38")
    
## ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "CD38.CD8.Naive", "CD38.CD8.CM", "CD38.CD8.EM", "CD38.CD8.EMRA"  ))
## subplot2(ggdf,"CD38")
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c("PD1.B","PD1.CD4","PD1.CD8"))
p7=subplot4(ggdf,"PD1","(C) percentage of PD1+ cells in sorted T cells" )
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "PD1.CD4.Naive", "PD1.CD4.CM", "PD1.CD4.EM", "PD1.CD4.EMRA"  ))
p8=subplot4(ggdf,"PD1","")
    
ggdf <- melt(dat0,id.vars = c("tissue","type", "Sex","new"),measure.vars=c( "PD1.CD8.Naive", "PD1.CD8.CM", "PD1.CD8.EM", "PD1.CD8.EMRA"  ))
p9=subplot4(ggdf,"PD1","")
    
multiplot(p1,p4,p7,p2,p5,p8,p3,p6,p9,cols=3)
dev.off()


orig.dat=dat

for(analyze.spleen in c(TRUE,FALSE))
  {

if(analyze.spleen)
  {
    dat=orig.dat[orig.dat[,"tissue"]=="spleen",]
    pdf(file="../results/flow_plot_SPL_problem_samples_removed.pdf")
  }else
{
  dat=orig.dat[orig.dat[,"tissue"]=="PBL",]
  pdf(file="../results/flow_plot_PBL_problem_samples_removed.pdf")
}



### generate tables of mean and sd of cell proportions for each tissue

## sumtab=NULL
## sumtab=cbind(sumtab,apply(dat[,6:ncol(dat)],2,mean,na.rm=T))
## sumtab=cbind(sumtab,apply(dat[,6:ncol(dat)],2,sd,na.rm=T))

## sumtab=cbind(sumtab,apply(dat[dat[,2]=="NZO",6:ncol(dat)],2,mean,na.rm=T))
## sumtab=cbind(sumtab,apply(dat[dat[,2]=="NZO",6:ncol(dat)],2,sd,na.rm=T))

## sumtab=cbind(sumtab,apply(dat[dat[,2]=="B6",6:ncol(dat)],2,mean,na.rm=T))
## sumtab=cbind(sumtab,apply(dat[dat[,2]=="B6",6:ncol(dat)],2,sd,na.rm=T))

## sumtab=cbind(sumtab,apply(dat[dat[,3]=="F",6:ncol(dat)],2,mean,na.rm=T))
## sumtab=cbind(sumtab,apply(dat[dat[,3]=="F",6:ncol(dat)],2,sd,na.rm=T))

## sumtab=cbind(sumtab,apply(dat[dat[,3]=="M",6:ncol(dat)],2,mean,na.rm=T))
## sumtab=cbind(sumtab,apply(dat[dat[,3]=="M",6:ncol(dat)],2,sd,na.rm=T))

## sumtab=round(sumtab,1)
## colnames(sumtab)=c("total mean","sd","NZO mean","sd","B6 mean","sd","F mean","sd","M mean","sd")
## if(analyze.spleen) write.csv(sumtab,file="../results/sumtab_SPL.csv") else write.csv(sumtab,file="../results/sumtab_PBL.csv")




if(TRANSFORM) for(i in 6:ncol(dat)) dat[,i]=transform(as.numeric(as.matrix(dat[,i])))

##### do regression of cell proportion w.r.t. age and generate a table of coef and p-values  #####

## percent.vs.age=matrix(NA,nr=ncol(dat)-6+2,nc=3)
## for(i in 6:ncol(dat))
##   {
##     percent.vs.age[i-5,]=c(change.column.name(colnames(dat)[i]),signif(summary(lm(dat[,i]~dat$age+dat$Sex+dat$type))$coef[2,c(1,4)],2))
##   }

## asterick=rep("",nrow(percent.vs.age))
## asterick[p.adjust(as.numeric(percent.vs.age[,3]),"fdr")<0.05]="*"
## percent.vs.age=cbind( asterick  ,percent.vs.age)
## colnames(percent.vs.age)=c("","celltype","coef","p")
## if(analyze.spleen) write.csv(percent.vs.age,file="../results/flow_percent_vs_age_spleen.csv",quote=F,row.names=F) else write.csv(percent.vs.age,file="../results/flow_percent_vs_age_PBL.csv",quote=F,row.names=F)

##### investigate effect of gender and strain on aging pattern #####
### test if the intercept and slope for aging effect is different between gender and strain ### For the column "M-F in B6” , 1 means slope/intercept for male is larger than female, -1 means the opposite and blank means that there is no significant difference

temp=summary(lm(dat[,6]~dat$age*dat$Sex*dat$type))$coef
percent.vs.type.interaction.coef=percent.vs.type.interaction.p=matrix(NA,nr=ncol(dat)-5,nc=nrow(temp))
colnames(percent.vs.type.interaction.coef)=colnames(percent.vs.type.interaction.p)=rownames(temp)
for(i in 6:ncol(dat))
  {
    temp=signif(summary(lm(dat[,i]~dat$age*dat$Sex*dat$type))$coef[,1],2)
    percent.vs.type.interaction.coef[i-5,names(temp)]=temp

    temp=signif(summary(lm(dat[,i]~dat$age*dat$Sex*dat$type))$coef[,4],2)
    percent.vs.type.interaction.p[i-5,names(temp)]=temp
        
  }
rownames(percent.vs.type.interaction.p)=colnames(dat)[6:ncol(dat)]


tempdat=as.matrix(dat)
tempdat[tempdat[,"type"]=="NZO","type"]="ANZO"
tempdat[tempdat[,"Sex"]=="M","Sex"]="AM"
tempdat2=dat
tempdat2$"type"=tempdat[,"type"]
tempdat2$"Sex"=tempdat[,"Sex"]
tempdat=tempdat2


temp=summary(lm(tempdat[,6]~tempdat$age*tempdat$Sex*tempdat$type))$coef
percent.vs.type.interaction.coef2=percent.vs.type.interaction.p2=matrix(NA,nr=ncol(dat)-5,nc=nrow(temp))
colnames(percent.vs.type.interaction.coef2)=colnames(percent.vs.type.interaction.p2)=rownames(temp)
for(i in 6:ncol(dat))
  {
    temp=signif(summary(lm(tempdat[,i]~tempdat$age*tempdat$Sex*tempdat$type))$coef[,1],2)
    percent.vs.type.interaction.coef2[i-5,names(temp)]=temp

    temp=signif(summary(lm(tempdat[,i]~tempdat$age*tempdat$Sex*tempdat$type))$coef[,4],2)
    percent.vs.type.interaction.p2[i-5,names(temp)]=temp        
  }
rownames(percent.vs.type.interaction.p2)=colnames(dat)[6:ncol(dat)]

                


### variables that have significant interaction 
 temp=rowSums(percent.vs.type.interaction.p[,c("dat$SexM","dat$typeNZO","dat$age:dat$SexM","dat$age:dat$typeNZO","dat$age:dat$SexM:dat$typeNZO")]<0.005)>0
 sel=percent.vs.type.interaction.p[temp,1]
# dat0=dat[,c( "tissue","type","Sex","age","Sample", sel )]
dat0=dat


    datF=dat0[dat0$Sex=="F",]
    datM=dat0[dat0$Sex=="M",]
   
    datB6=dat0[dat0$type=="B6",]
    datNZO=dat0[dat0$type=="NZO",]
   

    datFB6=dat0[dat0$Sex=="F"& dat0$type=="B6",]
    datMNZO=dat0[dat0$Sex=="M"&dat0$type=="NZO",]

    datFNZO=dat0[dat0$Sex=="F"& dat0$type=="NZO",]
    datMB6=dat0[dat0$Sex=="M"&dat0$type=="B6",]


    plotting2(dat0,datMB6,datMNZO,datFB6,datFNZO)
    

#### fits regression

R1=regression(datFB6)
R2=regression(datMB6)
R3=regression(datFNZO)
R4=regression(datMNZO)

R5=regression(datB6)
R6=regression(datNZO)
if(analyze.spleen)
  {
    B6.spleen=R5
    NZO.spleen=R6
  }else{
    B6.PBL=R5
    NZO.PBL=R6
  }


percent.vs.type.interaction.p= cbind(percent.vs.type.interaction.p2,percent.vs.type.interaction.p)
percent.vs.type.interaction.coef= cbind(percent.vs.type.interaction.coef2,percent.vs.type.interaction.coef)


### test if the intercept and slope for aging effect is different between gender and strain. For the column "M-F in B6” , 1 means slope/intercept for male is larger than female, -1 means the opposite and blank means that there is no significant difference

cd38var=c( "CD38.B", "CD38.CD4", "CD38.CD8", "CD38.CD4.Naive", "CD38.CD4.CM" , "CD38.CD4.EM" , "CD38.CD4.EMRA", "CD38.CD8.Naive", "CD38.CD8.CM" , "CD38.CD8.EM" ,  "CD38.CD8.EMRA" )

var=c("dat$SexM","tempdat$SexF","dat$typeNZO","tempdat$typeB6")
x=sign(percent.vs.type.interaction.coef[,var])*percent.vs.type.interaction.p[,var]
x[,c(2,4)]= -x[,c(2,4)]
colnames(x)=c("M-F in B6","M-F in NZO","NZO-B6 in F","NZO-B6 in M")
rownames(x)=rownames(R1)
#if(!analyze.spleen) x[cd38var,-1]=1

x=x[c(4:6,53,7:47),]
MAX=max(abs(-sign(x)*log(abs(x))),na.rm=T);breaksList = seq(-MAX, MAX, by = 1)
pheatmap(-sign(x)*log(abs(x)),cluster_cols = FALSE,cluster_rows = FALSE,scale="none",main="Difference of intercepts",color = colorRampPalette(c("blue","white","red"))(length(breaksList)),breaks = breaksList)
y=sign(x)*pvalue.convert(abs(x))
y[y==0]=""
if(analyze.spleen)
  {
    write.csv(y,file="../results/difference_of_intercept_SPL.csv",quote=F)
  } else
{
  write.csv(y,file="../results/difference_of_intercept_PBL.csv",quote=F)
}



maineffect.coef=cbind(R1[,"aging_coef"],R2[,"aging_coef"],R3[,"aging_coef"],R4[,"aging_coef"])
maineffect.p=cbind(R1[,"aging_p"],R2[,"aging_p"],R3[,"aging_p"],R4[,"aging_p"])
x=sign(maineffect.coef)*maineffect.p
colnames(x)=c("F in B6","M in B6","F in NZO","M in NZO")

coef.heatmap(x[c(4:6,53,7:47),],"Slope with aging")

var=c("dat$age:dat$SexM","tempdat$age:tempdat$SexF" ,"dat$age:dat$typeNZO","tempdat$age:tempdat$typeB6")

x=sign(percent.vs.type.interaction.coef[,var])*percent.vs.type.interaction.p[,var]
x[,c(2,4)]= -x[,c(2,4)]
colnames(x)=c("M-F in B6","M-F in NZO","NZO-B6 in F","NZO-B6 in M")
rownames(x)=rownames(R1)
#if(!analyze.spleen) x[cd38var,-1]=1

coef.heatmap(x[c(4:6,53,7:47),],"Difference of slopes")
y=sign(x)*pvalue.convert(abs(x))
y[y==0]=""


if(analyze.spleen) write.csv(y,file="../results/difference_of_slope_SPL.csv",quote=F) else write.csv(y,file="../results/difference_of_slope_PBL.csv",quote=F)

dev.off()

}


ggpoint <- function(x,y,xl,yl,title)
  {
dat=data.frame(x=x,y=y,col=names(x))
p <- ggplot(dat, aes(x,y,color=col))+ geom_point()+labs(x=xl,y=yl,title=title)
return(p)
}
pdf(file="../results/FC_comparison_flow.pdf",width=10)
sel=4:14
x1=ggpoint(B6.spleen[sel,"aging_coef"],NZO.spleen[sel,"aging_coef"],"B6","NZO","spleen")
x2=ggpoint(B6.PBL[sel,"aging_coef"],NZO.PBL[sel,"aging_coef"],"B6","NZO","PBL")
x3=ggpoint(B6.spleen[sel,"aging_coef"],B6.PBL[sel,"aging_coef"],"spleen","PBL","B6")
x4=ggpoint(NZO.spleen[sel,"aging_coef"],NZO.PBL[sel,"aging_coef"],"spleen","PBL","NZO")
multiplot(x1,x2,x3,x4,cols=2)
sel=15:47
dev.off()

