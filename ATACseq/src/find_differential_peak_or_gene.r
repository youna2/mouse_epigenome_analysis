
p.cutoff=0.05

p.mat=NULL ##bed[,"p.value"]
fc.mat= NULL  ##-bed[,"Fold"]

#### remove samples of age 26.75 months ######

if(selB6) selectsample= (TYPE=="B6" & AGE != 26.75) else selectsample= (TYPE=="NZO" & AGE != 26.75)
y0=Y[,selectsample]
y0forheatmap=bed[,selectsample]

age0=AGE[selectsample]
gender0=GENDER[selectsample]
tissue0=TISSUE[selectsample]
type0=TYPE[selectsample]
samplemouseid0=SAMPLEMOUSEID[selectsample]

## Do timeseries analysis separately for each tissue type and gender #######
setwd("../results")
utissue=c( "spleen", "BM", "memory","naive", "PBL")         

if(selB6)
  {
    pdf(file="diffpeakB6.pdf")
    topgene="B6topgene"
  }else{
    pdf(file="diffpeakNZO.pdf")
    topgene="NZOtopgene"
  }
twoway.barplot.argF1=twoway.barplot.argF2=twoway.barplot.argF3=NULL
twoway.barplot.argM1=twoway.barplot.argM2=twoway.barplot.argM3=NULL
for(i in 1:length(utissue))
  {
    y=y0[,tissue0==utissue[i]& gender0=="M"]
    age=age0[tissue0==utissue[i]& gender0=="M"]
    gender=gender0[tissue0==utissue[i]& gender0=="M"]
    type=type0[tissue0==utissue[i]& gender0=="M"]

### Do timeseries analysis within male
    atac.glmtopM=edgeRfit(y,age,gender,type)
    
    
    y=y0[,tissue0==utissue[i]& gender0=="F"]
    age=age0[tissue0==utissue[i]& gender0=="F"]
    gender=gender0[tissue0==utissue[i]& gender0=="F"]
    type=type0[tissue0==utissue[i]& gender0=="F"]

### Do timeseries analysis within female
    atac.glmtopF=edgeRfit(y,age,gender,type)
    
#### print number of age-increasing/decreasing peaks or genes  #################
    print("In Tissue") 
    print(utissue[i])

    print("In Female, significantly changing peaks/genes")
    print( sum(atac.glmtopF[,"FDR"]<p.cutoff))
    print("In Female, significantly increasing peaks/genes")
    print( sum(atac.glmtopF[,"FDR"]<p.cutoff & atac.glmtopF[,"logFC"]>0))
    print("In Female, significantly decreasing peaks/genes")
    print( sum(atac.glmtopF[,"FDR"]<p.cutoff & atac.glmtopF[,"logFC"]<0))

    twoway.barplot.argF1=c(twoway.barplot.argF1,rep(paste(utissue[i]),2))
    twoway.barplot.argF2=c(twoway.barplot.argF2,c( sum(atac.glmtopF[,"FDR"]<p.cutoff & atac.glmtopF[,"logFC"]>0), -sum(atac.glmtopF[,"FDR"]<p.cutoff & atac.glmtopF[,"logFC"]<0)))
    twoway.barplot.argF3=c(twoway.barplot.argF3,c("+","-"))
      
    
    print("In Male, significantly changing peaks/genes")
    print( sum(atac.glmtopM[,"FDR"]<p.cutoff))
    print("In Male, significantly increasing peaks/genes")
    print( sum(atac.glmtopM[,"FDR"]<p.cutoff & atac.glmtopM[,"logFC"]>0))
    print("In Male, significantly decreasing peaks/genes")
    print( sum(atac.glmtopM[,"FDR"]<p.cutoff & atac.glmtopM[,"logFC"]<0))

    twoway.barplot.argM1=c(twoway.barplot.argM1,rep(paste(utissue[i]),2))
    twoway.barplot.argM2=c(twoway.barplot.argM2,c( sum(atac.glmtopM[,"FDR"]<p.cutoff & atac.glmtopM[,"logFC"]>0), -sum(atac.glmtopM[,"FDR"]<p.cutoff & atac.glmtopM[,"logFC"]<0)))
    twoway.barplot.argM3=c(twoway.barplot.argM3,c("+","-"))
      
    
    print("peaks/genes that are significantly changing in both male and female")
    print( sum(atac.glmtopM[,"FDR"]<p.cutoff & atac.glmtopF[,"FDR"]<p.cutoff))

#### draw heatmap of age-increasing/decreasing peaks or genes  #################

    for(sex in c("F","M"))
      {
        if(sex=="M") atac.glmtop=atac.glmtopM
        if(sex=="F") atac.glmtop=atac.glmtopF

        p.mat=cbind(p.mat,atac.glmtop[,"FDR"])
        fc.mat=cbind(fc.mat,atac.glmtop[,"logFC"])

    if(sum(atac.glmtop[,"FDR"]<p.cutoff)>2)
      {

        heatmapmat=y0forheatmap[atac.glmtop[,"FDR"]<p.cutoff,tissue0==utissue[i]& gender0==sex]

        annot.row=annot=atac.glmtop[atac.glmtop[,"FDR"]<p.cutoff,"logFC"]
        annot[annot.row>0]="opening"
        annot[annot.row<0]="closing"
    
        colnames(heatmapmat)=paste(  samplemouseid0[tissue0==utissue[i]& gender0==sex] ,"_",age0[tissue0==utissue[i]  &  gender0==sex ],sep="")
        heatmapmat=heatmapmat[,order(age0[tissue0==utissue[i] & gender0==sex])]

        if(nrow(heatmapmat)>40000)
          {
            tmpsample=sample(1:nrow(heatmapmat),40000)
            heatmapmat=heatmapmat[tmpsample,]
            annot=annot[tmpsample]
          }

        annot=as.data.frame(annot)
        rownames(annot)=rownames(heatmapmat)

        if(min(heatmapmat)==0) heatmapmat=heatmapmat+1



        pheatmap(log(heatmapmat),scale="row",cluster_cols = FALSE,main=paste(utissue[i],sex, sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]>0),"opening", sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]<0),"closing"),annotation_row=annot,show_rownames=F,color=colorRampPalette(c("blue","white","red"))(100))
      }
        
      }

  }
colnames(fc.mat)=colnames(p.mat)=paste(rep(utissue,each=2),rep(c("F","M"),length(utissue)))
save(p.mat,fc.mat,file=paste("pmat_fcmat_B6_",selB6,".Rdata",sep=""))
### heatmap of p-values of peaks/genes across tissues and gender
global.heatmap(p.mat,fc.mat)



### peaks/genes that are commonly increasing/decreasing across tissues and gender
common.peaks(p.mat,fc.mat,TRUE,topgene,annotation)

f.increasing=nrow(read.delim(paste(topgene,"_F_increasing",".txt",sep="")))
f.decreasing=nrow(read.delim(paste(topgene,"_F_decreasing",".txt",sep="")))
m.increasing=nrow(read.delim(paste(topgene,"_M_increasing",".txt",sep="")))
m.decreasing=nrow(read.delim(paste(topgene,"_M_decreasing",".txt",sep="")))

twoway.barplot.argF1=c(twoway.barplot.argF1,rep("common",2))
twoway.barplot.argF2=c(twoway.barplot.argF2,c(f.increasing,-f.decreasing))
twoway.barplot.argF3=c(twoway.barplot.argF3,c("+","-"))

twoway.barplot.argM1=c(twoway.barplot.argM1,rep("common",2))
twoway.barplot.argM2=c(twoway.barplot.argM2,c(m.increasing,-m.decreasing))
twoway.barplot.argM3=c(twoway.barplot.argM3,c("+","-"))

ylimmax=max(abs(c(twoway.barplot.argM2,twoway.barplot.argF2)))
ylimmax=40000;YLIM=c(-ylimmax,ylimmax)

q1=twoway.barplot(twoway.barplot.argF1,twoway.barplot.argF2,twoway.barplot.argF3,(-10):10*1000,(-10):10*1000,"Tissue","no. differential peaks/genes",paste(type0[1]," F",sep=""),YLIM)

q2=twoway.barplot(twoway.barplot.argM1,twoway.barplot.argM2,twoway.barplot.argM3,(-10):10*1000,(-10):10*1000,"Tissue","no. differential peaks/genes",paste(type0[1]," M",sep=""),YLIM)

multiplot(q1,NA,q2,NA,cols=2)



tissue.gender.type <- c(colnames(p.mat),"_F_increasing","_F_decreasing","_M_increasing","_M_decreasing","_all_increasing","_all_decreasing")

### save differential peaks/genes of each tissue as a txt file ####
diff.peaks(p.mat,fc.mat,topgene)

### See if the age-related pattern is common across tissues

    for(jj in 1:2)
      {
        if(jj==1) tt=((p.mat<p.cutoff & fc.mat>0)) else tt=((p.mat<p.cutoff & fc.mat<0))

        fisher.p=fisher.p0=fisher.stat=matrix(NA,nr=ncol(p.mat),nc=ncol(p.mat))
        rownames(fisher.p)=colnames(fisher.p)=colnames(p.mat)
        rownames(fisher.stat)=colnames(fisher.stat)=colnames(p.mat)

        for(i in 1:ncol(tt))
          for(j in 1:ncol(tt))
            {
              print(c(colnames(tt)[i],colnames(tt)[j]));
              print(table(tt[,i],tt[,j]));
              temp=table(tt[,i],tt[,j])
              if(ncol(temp)>1 & nrow(temp)>1)
                {
                  total=sum(temp)
                  black=sum(temp[1,])
                  white=sum(temp[2,])
                  pick=sum(temp[,2])
                  whitepick=temp[2,2]-1
                  fisher.p0[i,j]=  1-phyper(whitepick,white,black,pick)
                  fisher.p[i,j]=   fisher.test(temp,alternative="greater")$"p.value"
                  fisher.stat[i,j]=   fisher.test(temp,alternative="greater")$"estimate"
                }
            }
        print(mean(abs(fisher.p-fisher.p0)),na.rm=T) ## to check if my calculation is correct
        fisher.p=signif(fisher.p,2)
        fisher.p[upper.tri(fisher.p,diag=T)]="*"
        fisher.stat[upper.tri(fisher.stat,diag=T)]= 0

        if(jj==1)
          {
            write.csv(fisher.p,file=paste("fisher_pvalue_increasing_B6_",selB6,".csv",sep=""),quote=F)

            pheatmap(fisher.stat,scale="none",cluster_cols = FALSE,cluster_rows=FALSE,main=paste("overlap of increasing peaks in",TYPE,"(odds ratio)"))
            #
          } else
        {
          write.csv(fisher.p,file=paste("fisher_pvalue_decreasing_B6_",selB6,".csv",sep=""),quote=F)
          
          pheatmap(fisher.stat,scale="none",cluster_cols = FALSE,cluster_rows=FALSE,main=paste("overlap of decreasing peaks in ",TYPE,"(odds ratio)"))
        }
      }
dev.off()  
#### Do pathway analysis using immune genes ####

library("biomaRt")

load("../../ATACseq/data/biomaRt_human_mouse.Rdata")

load("../../ATACseq/data/mousehumangene_annotation.Rdata")
all.path.res=vector("list",2)
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

enrichpath=vector("list",3)
for(i in 1:3) enrichpath[[i]]=vector("list",length(tissue.gender.type))
      
for(N in 1:3)
  {
    directionsel=N
for(k in 1:length(tissue.gender.type))
  {
    temptissue=tissue.gender.type[k]

    scanfile=scan(paste(topgene,temptissue,".txt",sep=""))

    if(length(scanfile)>2)
      {
   
        diff.gene=as.matrix(read.table(paste(topgene,temptissue,".txt",sep=""),header=F))### differential gene
        if(directionsel==2) diff.gene=rbind(diff.gene[diff.gene[,3]>0,])
        if(directionsel==3) diff.gene=rbind(diff.gene[diff.gene[,3]<0,])
        if(nrow(diff.gene)>1)
          {
        genesV2 = getLDS(attributes = c("entrezgene"), filters = "entrezgene", values = diff.gene[,1] , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      }
        if(nrow(genesV2)>1)
          {

        genesV2[,2]=toupper(genesV2[, 2])
        
        gene.and.CA=diff.gene[match(genesV2[,1],diff.gene[,1]),]

        gene.and.CA[,1]==genesV2[,1]


        human.diff.gene <- unique(genesV2[, 2]) ### human ortholog of the differential gene
        write.table(human.diff.gene,file=paste(topgene,temptissue,"human.txt",sep=""),quote=F,row.names=F,col.names=F)
        allpath=NULL

        mean(human.diff.gene %in% gene.universe)

          if(immunemodule)
            {
              version="2008"#"2015"
            all.gene=as.matrix(read.table(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_genes.txt",sep=""),header=T))
            path.annotation=as.matrix(read.csv(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_annotations.csv",sep="")))
          }
        if(celltype.annotation)
          {
            load("../../ATACseq/data/pbmc_specific_genes.annotations_April.2018.EJM_10x.RData")
            all.gene=  geneset.genes.scrnaseq_pbmc_specific
            path.annotation= geneset.names.scrnaseq_pbmc_specific[,c(2,1)]
          }
            pathid=unique(all.gene[,1])
            pathp=rep(NA,length(pathid))
            for(i in 1:length(pathid))
              {
                path.gene=toupper(all.gene[all.gene[,1]==pathid[i],2])### all genes in pathway i

                total=length(unique(gene.universe))  ### all human ortholog genes 
                white=length(unique(intersect(path.gene,gene.universe))) ### all genes in pathway i in universe
                black=total-white
                pick=length(human.diff.gene)

                intersection=intersect(human.diff.gene,path.gene)
                
                lintersection= length(intersection )

                whitepick=lintersection-1

                pathp[i]= 1-phyper(whitepick,white,black,pick)

                temp=match(intersection,genesV2[,2])
                
              }

            allpath=rbind(allpath,cbind(path.annotation[ match(pathid,path.annotation[,1]), 2],pathp))
          

        pick=rep(NA,nrow(allpath))
        for(i in 1:length(pick))
          pick[i]=pick_null_pathway(allpath[i,1])
        allpath=allpath[!pick,]
        print(paste("all pathways that are named with fdr<0.05 for",temptissue))
        enrichpath[[N]][[k]]=rbind(allpath[p.adjust(as.numeric(allpath[,2]),"fdr")<0.05 & allpath[,1]!="Unknown",])

        print(enrichpath[[N]][[k]])
      }
  }
  }
  }
    all.path.res[[pathwaytype]]=enrichpath
  }


#### draw plots of pathway analysis results ####

if(selB6)
  {
    save(all.path.res,file="enrichpathwayB6.Rdata")
    pdf(file="pathwayplotB6.pdf")
  }else{
    save(all.path.res,file="enrichpathwayNZO.Rdata")
    pdf(file="pathwayplotNZO.pdf")
  }

source("../../ATACseq/src/plot.r")
dev.off()
