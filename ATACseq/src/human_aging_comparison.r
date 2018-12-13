library(corrplot)
source("../../ATACseq/src/library_function.r")
#balloonplot=vector("list",4)
textsize=6
RNA=FALSE
firstcohort=TRUE
onlypromoter=TRUE
genderseparate=TRUE
change.memory.name <- function(x)
  {
    for(i in 1:length(x))
      {
        x[[i]][x[[i]]=="memory"]="CD8.memory"
        x[[i]][x[[i]]=="naive"]="CD8.naive"
      }
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
    multiplot(temp[[1]],NA,temp[[2]],NA,temp[[3]],NA,temp[[4]],NA,temp[[5]],NA,cols=ncol)
  }


ATACprocess <- function(x,title)
  {
                                        #    if(onlypromoter)  x=x[x[,"Annotation"]%in%"promoter-TSS",]
    if(onlypromoter)  x=x[abs(as.numeric(x[,"DistancetoTSS"]))<=1000,]
    
    x=extract.closest.peak(x)

    Contrast=rep(title,nrow(x))
    x=cbind(x,Contrast)   
  return(x)
}

if(RNA)
  {
    load( "/Users/youna/mouse_aging_data_analysis_code/RNAseq/results/balloonplotRNA.Rdata")
    balloonplot=change.memory.name(balloonplot)
    
    tid=100


    if(!firstcohort)
      {    
        tissue.gender.type=c("male human", "female human")
        tissue.gender.type2=c("Males_Age3xAge1", "Females_Age3xAge1")
        tissue.gender.type3=c("PBL", "PBL")

        genename=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/aging.summer2017_rna_expressed_plain.txt"))

        if(genderseparate) RNAresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/aging.summer2017_rna_glm.results_AgexSex_homo_batchdate_BySex.txt")) else RNAresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/secondcohort/aging.summer2017_rna_glm.results_AgexSex_homo_batchdate.txt"))

        RNAresult[,"EnsemblID"]=genename[match(RNAresult[,"EnsemblID"],genename[,"EnsemblID"]),"GeneName"]
        colnames(RNAresult)[colnames(RNAresult)=="EnsemblID"]="GeneName"
      }else{
        tissue.gender.type=c("PBL human")
        tissue.gender.type2=tissue.gender.type3=c("PBL")

        RNAresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/firstcohort/pbmc_rna_glm.results.txt"))        
      }


    Contrast=rep("PBL",nrow(RNAresult))
    RNAresult=cbind(RNAresult,Contrast)    
    humangene=RNAresult[,c("GeneName","logFC","PValue","FDR","Contrast")]
   
    
    load("../../RNAseq/results/mouseRNAFC.Rdata")
    pdf(file="summary_pathwayplot_RNA.pdf",height=17,width=10)
    

  }else{

    load( "/Users/youna/mouse_aging_data_analysis_code/ATACseq/results/balloonplotATAC.Rdata")
    balloonplot=change.memory.name(balloonplot)
    

    tid=101

    if(!firstcohort)
      {
        tissue.gender.type=c("male human", "female human")#,"memory human","naive human")
        tissue.gender.type2=c("Males_Age3xAge1", "Females_Age3xAge1")#,"memory","naive")
        tissue.gender.type3=c("PBL", "PBL")#,"memory","naive")

        peaktogene=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/aging.summer2017_narrow_whitelisted_filtered_annotated_plain.txt"))

        peaktogene=extract.closest.peak(peaktogene)
        
        if(genderseparate)
          {
            pblresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/aging.summer2017_narrowPeaks_glm.results_AgexSex_homo_batchdate_BySex.txt"))
          }else {
            pblresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/secondcohort/da_output_AgexSex_homo_batchdate.glm.results.txt"))
            Contrast=rep("PBL",nrow(pblresult))
            pblresult=cbind(pblresult,Contrast)
          }
        pblresult[,"peakco"]= peaktogene[match(pblresult[,"peakco"],peaktogene[,"peakco"]),"GeneName"]
        colnames(pblresult)[colnames(pblresult)=="peakco"]="GeneName"
        
      }else{
        tissue.gender.type=c("PBL human")#,"memory human","naive human")
        tissue.gender.type2=tissue.gender.type3=c("PBL")#,"memory","naive")

        pblresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm.results.txt"))
        load("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/firstcohort/unique_matching_peaks.Rdata")
        #pblresult <- ATACprocess(pblresult,"PBL")
        pblresult=pblresult[match(unique.matching.peaks,pblresult[,"peakco"]),]
        print(mean(unique.matching.peaks==pblresult[,"peakco"]))
        Contrast=rep("PBL",nrow(pblresult))
        pblresult=cbind(pblresult,Contrast)   
      }
    
    humangene=pblresult[,c("GeneName","logFC","PValue","FDR","Contrast")]
    
    ## cd8result=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/sortedcells/cd8_whitelisted_filtered_glm.results.txt"))
    ## cd8result <- ATACprocess(cd8result,"CD8")                
    ## humangene=rbind(humangene,cd8result[,c("GeneName","logFC","PValue","FDR","Contrast")])
    
    ## cd8memresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/sortedcells/cd8mem_whitelisted_filtered_glm.results.txt"))
    ## cd8memresult <- ATACprocess(cd8memresult,"CD8.memory")
    ## humangene=rbind(humangene,cd8memresult[,c("GeneName","logFC","PValue","FDR","Contrast")])
    
    ## cd8nairesult=ATACsortedresult=as.matrix(read.delim("~/mouse_aging_epigenome_research/ATACseq_data_analysis/aging_analysis_human_Eladio/sortedcells/cd8nai_whitelisted_filtered_glm.results.txt"))
    ## cd8nairesult <- ATACprocess(cd8nairesult,"CD8.naive")
    ## humangene=rbind(humangene,cd8nairesult[,c("GeneName","logFC","PValue","FDR","Contrast")])

    humangene=humangene[!is.na(humangene[,1]),]
    
 if(onlypromoter) load("../../ATACseq/results/mouseATACFCpromoter.Rdata") else load("../../ATACseq/results/mouseATACFC.Rdata")
    pdf(file="summary_pathwayplot_ATAC.pdf",height=17,width=10)
    
  }

NZOFC=fc.mat.NZO
B6FC=fc.mat.B6
colnames(NZOFC)=colnames(B6FC)=c("spleen","BM","CD8.memory","CD8.naive","PBL")
    

total.FC.graph=vector("list",2*length(tissue.gender.type))   
count=0
for(k in 1:length(tissue.gender.type2))
  {
    temp=humangene[humangene[,"Contrast"]%in% tissue.gender.type2[k] ,]
    
    commongenename=intersect(toupper(humangenename),toupper(temp[,1]))
    count=count+1
    total.FC.graph[[count]]=plot.correlationplot(B6FC[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,1])),"logFC"]),paste("B6 vs. human",tissue.gender.type2[k]),"log FC in B6","log FC in human")
    count=count+1
    total.FC.graph[[count]]= plot.correlationplot(NZOFC[match(commongenename,toupper(humangenename)),tissue.gender.type3[k]],as.numeric(temp[match(commongenename,toupper(temp[,1])),"logFC"]),paste("NZO vs. human",tissue.gender.type2[k]),"log FC in NZO","log FC in human")

  }
 multiplot(total.FC.graph[[1]],total.FC.graph[[2]],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,cols=3)# else multiplot(total.FC.graph[[1]],total.FC.graph[[2]],NA,NA,total.FC.graph[[3]],total.FC.graph[[4]],NA,NA,total.FC.graph[[5]],total.FC.graph[[6]],NA,NA,cols=3)



library("biomaRt")

load("../../ATACseq/data/biomaRt_human_mouse.Rdata")

load("../../ATACseq/data/mousehumangene_annotation.Rdata")
gene.universe=toupper(unique(humangene[,1]))

#gene.universe=human.gene.universe
setwd("../results")
previous.dir=getwd()

fdr.cutoff=0.1
all.path.res=vector("list",4)
for(pathwaytype in 1:2)
  {
    if(pathwaytype==1)
      {
        immunemodule=FALSE
        celltype.annotation=TRUE
      }else{
        immunemodule=TRUE
        celltype.annotation=FALSE
      }

    enrichpath=enrichpath.wiki=enrichpath.kegg=vector("list",3)
    for(i in 1:3) enrichpath[[i]]=enrichpath.wiki[[i]]=enrichpath.kegg[[i]]=vector("list",length(tissue.gender.type))
    
    for(N in 2:3)
      {
        directionsel=N
        for(k in 1:length(tissue.gender.type))
          {
            humangene2=humangene[humangene[,"Contrast"]%in% tissue.gender.type2[k] ,]

            if(directionsel==2)
              {
                human.diff.gene=humangene2[as.numeric(humangene2[,"logFC"])>0 & as.numeric(humangene2[,"FDR"])<fdr.cutoff,1]
                directionselsymbol="positive"
              }
                
            if(directionsel==3)
              {

                human.diff.gene=humangene2[as.numeric(humangene2[,"logFC"])<0 & as.numeric(humangene2[,"FDR"])<fdr.cutoff,1]
                directionselsymbol="negative"
              }

            source("../../ATACseq/src/pathway_enrichment_test.r")
              
          }
      }
  
all.path.res[[pathwaytype]]=enrichpath
}
all.path.res[[3]]=enrichpath.wiki
all.path.res[[4]]=enrichpath.kegg
source("../../ATACseq/src/balloonplot.r")
                                        # save(all.path.res,file="enrichpathwayStrain.Rdata")


plot.balloonplot.human(balloonplot)
plot.balloonplot.mouse(balloonplot)


if(RNA)
  {
    load("../../RNAseq/results/mouseRNAFC.Rdata")
  } else {
   load("../../ATACseq/results/mouseATACFC.Rdata")
 }

NZOFC=fc.mat.NZO
B6FC=fc.mat.B6
colnames(NZOFC)=colnames(B6FC)=c("spleen","BM","CD8.memory","CD8.naive","PBL")
    

par(mfrow=c(2,2))
corrplot(corFC.bet.celltype(B6FC),type="upper",diag=FALSE,tl.srt=45)
mtext("correlation of FC between tissues in B6", at=3, line=-8, cex=1)
corrplot(corFC.bet.celltype(NZOFC),type="upper",diag=FALSE,tl.srt=45)
mtext("correlation of FC between tissues in NZO", at=3, line=-8, cex=1)


dev.off()



if(RNA)
  {
    pdf(file="FC_comparison_RNA.pdf",width=17)
    load("../../RNAseq/results/mouseRNAFC.Rdata")
  } else {
   pdf(file="FC_comparison_ATAC.pdf",width=17)
   load("../../ATACseq/results/mouseATACFC.Rdata")### This includes only promoter peaks
 }

plotFC.bet.celltype(B6FC,"in B6",5)
plotFC.bet.celltype(NZOFC,"in NZO",5)
plotFC.bet.strain(B6FC,NZOFC,5)
dev.off()
