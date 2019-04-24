library(corrplot)
source("../../ATACseq/src/library_function.r")

using.bnmapper=FALSE
### generating mouseRNAFC.Rdata mouseATACFC.Rdata ####
source("../../ATACseq/src/fcmat.r")
newhumandata=F


#balloonplot=vector("list",4)
textsize=6
RNA=TRUE
firstcohort=TRUE
readCD8=FALSE
genderseparate=TRUE
find.expressed.gene <- function(RNAresult)
  {
    sel=RNAresult[,"Contrast"]=="Females_Age3xAge1" |RNAresult[,"Contrast"]=="Males_Age3xAge1"

    RNAresult=RNAresult[sel,]
    gene=unique(RNAresult[,"GeneName"])

    y=rep(NA,length(gene))
    for(i in 1:length(gene))
      {
        y[i]=mean(as.numeric(RNAresult[RNAresult[,"GeneName"]%in% gene[i],"logCPM"]))
      }

    unexpressedgene=gene[y<=0]
    return(unexpressedgene)
  }

change.memory.name <- function(x)
  {
    for(i in 1:length(x))
      {
        x[[i]][x[[i]]=="memory"]="CD8.memory"
        x[[i]][x[[i]]=="naive"]="CD8.naive"
      }
    return(x)
  }
convert.col.name <- function(x)
  {
    for(i in 1:length(x)) x[i]=strsplit(x[i]," ")[[1]][1]
    return(x)
  }

extract.closest.peak <- function(peaktogene)
  {
    uniquegene=unique(peaktogene[,"GeneName"])
    newpeaktogene=matrix(NA,nr=length(uniquegene),nc=ncol(peaktogene))
    colnames(newpeaktogene)=colnames(peaktogene)
    for(j in 1:length(uniquegene))
      {
        temp=rbind(peaktogene[peaktogene[,"GeneName"]%in%uniquegene[j],])
        temp.distance=abs(as.numeric(temp[,"DistancetoTSS"]))
        newpeaktogene[j,]=temp[which(temp.distance==min(temp.distance))[1],]
      }

    
    return(newpeaktogene)
  }

plotFC.bet.celltype <- function(x,title,ncol)
  {
    count=0
    temp=vector("list",10)
    for(i in 1:(ncol(x)-1))
      for(j in (i+1):ncol(x))
        {
          count=count+1
          temp[[count]]=  plot.correlationplot(x[,i],x[,j],paste(colnames(x)[i],"vs",colnames(x)[j],title),paste("log FC in",colnames(x)[i]),paste("log FC in",colnames(x)[j]))
        }
    multiplot(temp[[1]],temp[[2]],temp[[3]],temp[[4]],temp[[5]],temp[[6]],temp[[7]],temp[[8]],temp[[9]],temp[[10]],cols=ncol)
  }

corFC.bet.celltype <- function(x)
  {
    M=matrix(1,nr=ncol(x),nc=ncol(x))
    colnames(M)=rownames(M)=colnames(x)

    for(i in 1:ncol(x))
      for(j in 1:ncol(x))
        {
          M[i,j]= cor.test(x[,i],x[,j])$estimate
        }
    return(M)    
  }

plotFC.bet.strain <- function(x,y,ncol)
  {
    temp=vector("list",ncol(x))
    for(i in 1:ncol(x))
      {
          temp[[i]]=  plot.correlationplot(x[,i],y[,i],colnames(x)[i],"log FC in B6","log FC in NZO")
        }
    multiplot(temp[[1]],NA,NA,NA,temp[[2]],NA,NA,NA,temp[[3]],NA,NA,NA,temp[[4]],NA,NA,NA,cols=ncol)
  }


ATACprocess <- function(x,title)
  {
                                        #    if(onlypromoter)  x=x[x[,"Annotation"]%in%"promoter-TSS",]
    x=x[abs(as.numeric(x[,"DistancetoTSS"]))<=1000,]
    
    x=extract.closest.peak(x)

    Contrast=rep(title,nrow(x))
    x=cbind(x,Contrast)   
  return(x)
}

if(RNA)
  {
    load( "/Users/youna/mouse_aging_data_analysis_code/RNAseq/results/balloonplot.Rdata")
    balloonplot=change.memory.name(balloonplot)
    
    tid=100


    if(!firstcohort)
      {    
        tissue.gender.type=c("male human", "female human")
        tissue.gender.type2=c("Males_Age3xAge1", "Females_Age3xAge1")
        tissue.gender.type3=c("PBL", "PBL")

        genename=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/aging.summer2017_rna_expressed_plain.txt"))

        if(genderseparate) RNAresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/aging.summer2017_rna_glm.results_AgexSex_homo_batchdate_BySex.txt")) else RNAresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/secondcohort/aging.summer2017_rna_glm.results_AgexSex_homo_batchdate.txt"))

        RNAresult[,"EnsemblID"]=genename[match(RNAresult[,"EnsemblID"],genename[,"EnsemblID"]),"GeneName"]
        colnames(RNAresult)[colnames(RNAresult)=="EnsemblID"]="GeneName"

        if(genderseparate)
          {
            unexpressedgene=find.expressed.gene(RNAresult)
            save(unexpressedgene,file="../results/unexpressedgene_second_cohort.Rdata")
          }
      }else{
        tissue.gender.type="PBL human"
        tissue.gender.type2=tissue.gender.type3="PBL"

        RNAresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_rna_glm.results.txt"))        
      }


    Contrast=rep("PBL",nrow(RNAresult))
    RNAresult=cbind(RNAresult,Contrast)    
    humangene=RNAresult[,c("GeneName","logFC","PValue","FDR","Contrast","logCPM")]
   
    
    load("../../RNAseq/results/mouseRNAFC.Rdata")

    
    pdf(file="summary_pathwayplot_RNA.pdf",height=21,width=10)
    

  }else{

    load( "/Users/youna/mouse_aging_data_analysis_code/ATACseq/results/balloonplot.Rdata")
    balloonplot=change.memory.name(balloonplot)
    

    tid=101

    if(!firstcohort)
      {
        tissue.gender.type=c("male human", "female human")#,"memory human","naive human")
        tissue.gender.type2=c("Males_Age3xAge1", "Females_Age3xAge1")#,"memory","naive")
        tissue.gender.type3=c("PBL", "PBL")#,"memory","naive")

        peaktogene=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/aging.summer2017_narrow_whitelisted_filtered_annotated_plain.txt"))

        peaktogene=extract.closest.peak(peaktogene)
        
        if(genderseparate)
          {
            pblresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/aging.summer2017_narrowPeaks_glm.results_AgexSex_homo_batchdate_BySex.txt"))
          }else {
            pblresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/secondcohort/da_output_AgexSex_homo_batchdate.glm.results.txt"))
            Contrast=rep("PBL",nrow(pblresult))
            pblresult=cbind(pblresult,Contrast)
          }
        pblresult[,"peakco"]= peaktogene[match(pblresult[,"peakco"],peaktogene[,"peakco"]),"GeneName"]
        colnames(pblresult)[colnames(pblresult)=="peakco"]="GeneName"
        
      }else{

        if(readCD8)
          {
            tissue.gender.type=c("PBL human","memory human","naive human")
            tissue.gender.type2=tissue.gender.type3=c("PBL","CD8.memory","CD8.naive")
          }else{
            tissue.gender.type=c("PBL human")
            tissue.gender.type2=tissue.gender.type3=c("PBL")
          }
        
        ## pblresult=read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm.results.txt") ### This is the original human ATAC file from Eladio
        ## load("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/unique_matching_peaks.Rdata")#### These are selected peaks from Eladio that were mapped to expressed genes in distance 1500(?) bp
        
        ## #pblresult <- ATACprocess(pblresult,"PBL")
        ## pblresult=pblresult[match(unique.matching.peaks,pblresult[,"peakco"]),]

        ## write.table(pblresult,file="../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm_filtered.txt",quote=F,row.names=F,col.names=T,sep="\t")


        pblresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm_filtered.txt")) ### This is the human ATAC file after filitering selected peaks from Eladio

        
        Contrast=rep("PBL",nrow(pblresult))
        pblresult=cbind(pblresult,Contrast)   
      }
    
    humangene=pblresult[,c("GeneName","logFC","PValue","FDR","Contrast","logCPM","peakco")]


    ####new human data ######
    if(newhumandata)
      {
    expressed.gene.first.cohort=humangene[,"GeneName"]

    expressed.gene.second.cohort=read.delim("../../ATACseq/data/aging_analysis_human_Eladio/secondcohort2/aging.summer2017_RNA_annotations.txt")

    expressed.gene.second.cohort=expressed.gene.second.cohort[,"GeneName"]
    
    load("../results7/ATACseqData2.Rdata")
    load("../results7/pmat_fcmat_within B6_before_filtering.Rdata")
    
    mean(fc.mat.strain==atac.glmtop[,"logFC"])
    meanexp=atac.glmtop[,"logCPM"]
    
    humangene=cbind(annotation[,"Gene.Name"],fc.mat.strain,NA,p.mat.strain,"PBL",meanexp,annotation[,"peakco"])
    colnames(humangene)=c("GeneName","logFC","PValue","FDR","Contrast","logCPM","peakco")

    humangene=humangene[annotation[,"Distance.to.TSS"]<1500 & rowSums(bed>1)>4,]
    humangene=humangene[toupper(humangene[,"GeneName"])%in% toupper(expressed.gene.second.cohort),]
  }
    ## cd8result=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/sortedcells/cd8_whitelisted_filtered_glm.results.txt"))
    ## cd8result <- ATACprocess(cd8result,"CD8")                
    ## humangene=rbind(humangene,cd8result[,c("GeneName","logFC","PValue","FDR","Contrast")])



    if(readCD8)
      {
        cd8memresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/sortedcells/cd8mem_whitelisted_filtered_glm.results.txt"))
        cd8memresult <- ATACprocess(cd8memresult,"CD8.memory")
        humangene=rbind(humangene,cd8memresult[,c("GeneName","logFC","PValue","FDR","Contrast","logCPM","peakco")])
    
        cd8nairesult=ATACsortedresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/sortedcells/cd8nai_whitelisted_filtered_glm.results.txt"))
        cd8nairesult <- ATACprocess(cd8nairesult,"CD8.naive")
        humangene=rbind(humangene,cd8nairesult[,c("GeneName","logFC","PValue","FDR","Contrast","logCPM","peakco")])
      }

    
    humangene=humangene[!is.na(humangene[,1]),]
    
    
    load("../../ATACseq/results/mouseATACFC.Rdata")
    pdf(file="summary_pathwayplot_ATAC.pdf",height=21,width=10)
    
  }

NZOFC=fc.mat.NZO
B6FC=fc.mat.B6
if(RNA) tmpcolname=convert.col.name(colnames(NZOFC)) else tmpcolname="PBL"
colnames(meanexp.B6)=colnames(meanexp.NZO)=colnames(NZOFC)=colnames(B6FC)=colnames(p.mat.B6)=colnames(p.mat.NZO)=tmpcolname


#### draw FC comparison plot between human and mouse

nrow(meanexp.B6)==nrow(B6FC)

total.FC.graph=vector("list",2*length(tissue.gender.type))   
count=0
for(k in 1:length(tissue.gender.type2))
  {
    temp=humangene[humangene[,"Contrast"]%in% tissue.gender.type2[k] ,]

    if(using.bnmapper) PEAKID="peakco" else PEAKID="GeneName"
    
    commongenename=intersect(toupper(humangenename),toupper(temp[,PEAKID]))
    
    count=count+1
    total.FC.graph[[count]]=plot.correlationplot(B6FC[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,PEAKID])),"logFC"]),paste("B6 vs. human",tissue.gender.type2[k]),"log FC in B6","log FC in human")
    count=count+1
    total.FC.graph[[count]]= plot.correlationplot(NZOFC[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,PEAKID])),"logFC"]),paste("NZO vs. human",tissue.gender.type2[k]),"log FC in NZO","log FC in human")

  }


humanandmouseFC=cbind(temp[match(commongenename,toupper(temp[,PEAKID])),c("GeneName","logFC")],B6FC[match(commongenename,toupper(humangenename)),"PBL"],NZOFC[match(commongenename,toupper(humangenename)),"PBL"])

common.decreasing= rowMeans(matrix(as.numeric(humanandmouseFC[,-1]),nr=nrow(humanandmouseFC))<0)==1
common.increasing= rowMeans(matrix(as.numeric(humanandmouseFC[,-1]),nr=nrow(humanandmouseFC))>0)==1


humanandmouseFDR=cbind(temp[match(commongenename,toupper(temp[,PEAKID])),c("GeneName","FDR")],p.mat.B6[match(commongenename,toupper(humangenename)),"PBL"],p.mat.NZO[match(commongenename,toupper(humangenename)),"PBL"])

common.significant= rowMeans(matrix(as.numeric(humanandmouseFDR[,-1]),nr=nrow(humanandmouseFDR))<0.1)==1
#common.significant= as.numeric(humanandmouseFDR[,2])<0.1
common.significant= as.numeric(humanandmouseFDR[,2])<0.1  &  (as.numeric(humanandmouseFDR[,3])<0.1 |as.numeric(humanandmouseFDR[,4])<0.1 )
sum(common.decreasing)
sum(common.increasing)


common.increasing= common.increasing & common.significant
common.decreasing= common.decreasing & common.significant

sum(common.decreasing)
sum(common.increasing)

if(RNA) foldername="RNAseq" else foldername="ATACseq"
write.table(humanandmouseFC[common.increasing,1],file=paste("../../",foldername,"/results/common_increasing_genes.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(humanandmouseFC[common.decreasing,1],file=paste("../../",foldername,"/results/common_decreasing_genes.txt",sep=""),quote=F,row.names=F,col.names=F)


if(RNA | !readCD8) multiplot(total.FC.graph[[1]],total.FC.graph[[2]],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,cols=3) else multiplot(total.FC.graph[[1]],total.FC.graph[[2]],NA,NA,total.FC.graph[[3]],total.FC.graph[[4]],NA,NA,total.FC.graph[[5]],total.FC.graph[[6]],NA,NA,cols=3)

    
#### draw mean comparison plot between human and mouse

total.graph=vector("list",2*length(tissue.gender.type))   
count=0
for(k in 1:length(tissue.gender.type2))
  {
    temp=humangene[humangene[,"Contrast"]%in% tissue.gender.type2[k] ,]
    
    commongenename=intersect(toupper(humangenename),toupper(temp[,PEAKID]))
    count=count+1
    total.graph[[count]]=plot.correlationplot(meanexp.B6[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,PEAKID])),"logCPM"]),paste("B6 vs. human",tissue.gender.type2[k]),"log CPM in B6","log CPM in human")
    count=count+1
    total.graph[[count]]= plot.correlationplot(meanexp.NZO[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,PEAKID])),"logCPM"]),paste("NZO vs. human",tissue.gender.type2[k]),"log CPM in NZO","log CPM in human")

  }
if(RNA | !readCD8) multiplot(total.graph[[1]],total.graph[[2]],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,cols=3) else multiplot(total.graph[[1]],total.graph[[2]],NA,NA,total.graph[[3]],total.graph[[4]],NA,NA,total.graph[[5]],total.graph[[6]],NA,NA,cols=3)




source("../../ATACseq/src/pathwayanalysis.r") ### This code does the pathway enrichment analysis for immune modules, and kegg and wiki



source("../../ATACseq/src/balloonplot.r")  ### This code draws the pathway enrichment plot

source("../../ATACseq/src/barplot.r")  ### This code draws the number corresponding to cell types plot



source("../../ATACseq/src/pathwayanalysis_common.r") ### This code does the pathway enrichment analysis for immune modules, and kegg and wiki

tissue.gender.type="PBL common"

source("../../ATACseq/src/balloonplot.r")  ### This code draws the pathway enrichment plot


source("../../ATACseq/src/barplot.r")  ### This code draws the number corresponding to cell types plot


                                        # save(all.path.res,file="enrichpathwayStrain.Rdata")
#plot.balloonplot.human(balloonplot)
plot.balloonplot.mouse(balloonplot)
                                      
#
pathway.enrich.mat=NULL
for(i in c(1,2,7))
{
tmp=NULL
tmp=c(tmp,which(balloonplot[[i]][,5]=="within B6"))
tmp=c(tmp, which(balloonplot[[i]][,5]=="within NZO"))
tmp=c(tmp, which(balloonplot[[i]][,5]=="human"))
tmp=c(tmp, which(balloonplot[[i]][,5]=="common"))

pathway.enrich.mat=rbind(pathway.enrich.mat,cbind(c(pathwaytitle[i],rep("",length(tmp)-1)),balloonplot[[i]][tmp,]))
pathway.enrich.mat=rbind(pathway.enrich.mat,rep("",ncol(pathway.enrich.mat)))
}
write.csv(pathway.enrich.mat,file=paste("pathwayenrichment_RNA",RNA,".csv",sep=""),quote=F,row.names=F)
#
draw.barplot.human(barplot)
draw.barplot.mouse(barplot)
if(RNA) diff.genes() else diff.peaks()

### draw FC comparison plot between cell types in B6 and NZO 
    
if(RNA)
  {
   load("../../RNAseq/results/mouseRNAFC.Rdata")
 
   NZOFC=fc.mat.NZO
   B6FC=fc.mat.B6
   tmpcolname=convert.col.name(colnames(NZOFC))
   colnames(meanexp.B6)=colnames(meanexp.NZO)=colnames(NZOFC)=colnames(B6FC)=colnames(p.mat.B6)=colnames(p.mat.NZO)=tmpcolname

   plotFC.bet.strain(B6FC,NZOFC,4)

source("../src/comparison_B6_NZO_cor.r")
#p1
                                        #+ scale_fill_manual(  values=c("PBL"="red","spleen"= "blue"))+ scale_x_discrete(limits = rev(levels(factor(type))))



   
    ## par(mfrow=c(2,2))
    ## corrplot(corFC.bet.celltype(B6FC),type="upper",diag=FALSE,tl.srt=45)
    ## mtext("correlation of FC between tissues in B6", at=3, line=-8, cex=1)
    ## corrplot(corFC.bet.celltype(NZOFC),type="upper",diag=FALSE,tl.srt=45)
    ## mtext("correlation of FC between tissues in NZO", at=3, line=-8, cex=1)


}  else {

    temp=temp2=vector("list",5)
cormat=NULL

   for(BTID in 2:5)
      {
        load(paste("../../ATACseq/results",BTID,"/pmat_fcmat_within B6_before_filtering.Rdata",sep=""))
        B6FC=fc.mat.strain


        load(paste("../../ATACseq/results",BTID,"/pmat_fcmat_within NZO_before_filtering.Rdata",sep=""))
        NZOFC=fc.mat.strain

        temp[[BTID]]=  plot.correlationplot(B6FC[,1],NZOFC[,1],convert.col.name(colnames(B6FC)[1]),"log FC in B6","log FC in NZO")
humangenename=annotation.orig[,"Gene.Name"]
source("../src/comparison_B6_NZO_cor.r")
        cormat=rbind(cormat, df)
        
      }
    
    multiplot(temp[[2]],NA,NA,NA,temp[[3]],NA,NA,NA,temp[[4]],NA,NA,NA,temp[[5]],NA,NA,NA,cols=4)
df=cormat
 #       ggplot(data=df, aes(x=module, y=cor)) +  geom_bar(stat="identity", position=position_dodge())+  theme_minimal()+ggtitle("B6 vs. NZO")+coord_flip()+labs(x="")+facet_wrap(~tissue,nrow=1)


  }


### draw FC comparison plot between cell types 

## if(RNA)
##   {
##     pdf(file="FC_comparison_RNA.pdf",width=17)
##   } else {
##    pdf(file="FC_comparison_ATAC.pdf",width=17)
##  }

## plotFC.bet.celltype(B6FC,"in B6",4)
## plotFC.bet.celltype(NZOFC,"in NZO",4)

## dev.off()



znormalize <- function(x)
  {
    for(i in 1:nrow(x))
      x[i,]=(as.numeric(x[i,])-mean(as.numeric(x[i,])) )/sd(as.numeric(x[i,]))
    return(x)
    
  }

remove.dash <- function(x)
  {
    for(i in 1:length(x)) x[i]=strsplit(x[i],"-")[[1]][1]
    return(x)
  }
#load("../../RNAseq/results/PBL_expression_B6.Rdata")

if(RNA)
  {
    load("../../RNAseq/results/RNAseqData2.Rdata")
    humangenename=convert.mouse.entrez.to.human.symbol(annotation.orig[,1])
  }else{
    load("../../ATACseq/results3/ATACseqData2.Rdata")
    humangenename=convert.mouse.entrez.to.human.symbol(annotation.orig[,"Entrez.ID"])
  }


load("../../ATACseq/results/list_of_common_genes.Rdata")

gene=RNAlist$negative$PBL.PBMC.common
#genelist= gene[gene[,"T Cells"]==1,1]
genelist=gene[,1]
gene=RNAlist$positive$PBL.PBMC.common
genelist=c(genelist,gene[,1])

p=vector("list",150)
count=0

## human
if(RNA)
  {
    RNAresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_rna_glm.results.txt"))        
    RNAresult[,1:38]=znormalize(RNAresult[,1:38])
    humandata=RNAresult
    oldcolumn=1:16
    youngcolumn=17:38
  }else{
    pblresult=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm_filtered.txt")) ### This is the human ATAC file after filitering se
    pblresult[,11:47]=znormalize(pblresult[,11:47])
    humandata=pblresult
    oldcolumn=11:27
    youngcolumn=28:47
  }
humandata[,"GeneName"]=remove.dash(humandata[,"GeneName"])
## mouse
bed[,TISSUE=="PBL" & TYPE=="B6"]=znormalize(bed[,TISSUE=="PBL" & TYPE=="B6"])
bed[,TISSUE=="PBL" & TYPE=="NZO"]=znormalize(bed[,TISSUE=="PBL" & TYPE=="NZO"])

for(selectedgene in genelist)
  {
    count=count+1
    sel=which(humangenename==selectedgene)[1]
    sel2=which(toupper(humandata[,"GeneName"])==selectedgene)[1]

    if(length(sel)>0 & length(sel2)>0)
      {
    tmp=rbind(cbind(humandata[sel2,oldcolumn],"O","human"),cbind(humandata[sel2,youngcolumn],"Y","human"),cbind(bed[sel,TISSUE=="PBL"],AGE[TISSUE=="PBL"],TYPE[TISSUE=="PBL"]))

    tmp[tmp[,2]==18,2]="O"
    tmp[tmp[,2]==3,2]="Y"
    df=data.frame(exp=as.numeric(tmp[,1]),age=tmp[,2],strain=tmp[,3])
    df$age=factor(df$age,levels=c("Y","O"),ordered=TRUE)
    df$strain=factor(df$strain,levels=c("human","B6","NZO"),ordered=TRUE)
    p[[count]]=ggplot(df, aes(x=strain, y=exp, fill=age)) + geom_boxplot()+ggtitle(selectedgene)+xlab("")+ylab("expression level (z score)")
  }
  }
for(k in 0:14) multiplot(p[[1+k*10]],p[[2+k*10]],p[[3+k*10]],p[[4+k*10]],p[[5+k*10]],p[[6+k*10]],p[[7+k*10]],p[[8+k*10]],p[[9+k*10]],p[[10+k*10]],cols=2)



dev.off()
