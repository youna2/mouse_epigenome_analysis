
p.cutoff=0.05

p.mat.strain=NULL ##bed[,"p.value"]
fc.mat.strain= NULL  ##-bed[,"Fold"]
if(tid>=5)
  {
    commonpattern=TRUE
    if(tid==5) int.or.slope="common"
    if(tid==6) int.or.slope="opposite"
    if(tid==7) int.or.slope="NZO.signif"
    if(tid==8) int.or.slope="B6.signif"
    if(tid==9) int.or.slope="withinB6"
    if(tid==10) int.or.slope="withinNZO"
    if(tid==11) int.or.slope="all"
      
  }else{
    commonpattern=FALSE
    if(tid==3)
      {
        typeinteraction=TRUE
        int.or.slope="slope"
      }else{
        typeinteraction=FALSE
        int.or.slope="intercept"
      }
  }
#### remove samples of age 26.75 months ######

selectsample= (AGE != 26.75)
y0=Y[,selectsample]
y0forheatmap=bed[,-(1:4)][,selectsample]
age0=AGE[selectsample]
gender0=GENDER[selectsample]
tissue0=TISSUE[selectsample]
type0=TYPE[selectsample]
samplemouseid0=SAMPLEMOUSEID[selectsample]

## Do differential analysis between NZO and B6 in each tissue type ######

setwd("../results")
previous.dir=getwd()
utissue=c( "spleen", "BM", "memory","naive", "PBL")         

if(commonpattern)
  {
    pdf(paste(int.or.slope,"pattern.pdf",sep="_"))
    topgene=paste("B6_NZO",int.or.slope,sep="_")
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

    positivepeak="higher in NZO"
    negativepeak="lower in NZO"

    if(typeinteraction)
      {
        pdf("diffpeakstrain_B6_vs_NZO_slope.pdf")
        topgene="B6_vs_NZO_slope"
      }else{
        pdf("diffpeakstrain_B6_vs_NZO_intercept.pdf")
        topgene="B6_vs_NZO_intercept"
      }
  }



twoway.barplot.arg1=twoway.barplot.arg2=twoway.barplot.arg3=NULL

for(i in 1:length(utissue))
  {
    y=y0[,tissue0==utissue[i]]
    age=age0[tissue0==utissue[i]]
    gender=gender0[tissue0==utissue[i]]
    type=type0[tissue0==utissue[i]]
    if(length(unique(age[type=="B6"]))>1 & length(unique(age[type=="NZO"]))>1)
      {
        ### Do differential or common analysis between NZO and B6

        
        if(commonpattern)  atac.glmtop.strain=edgeRfitcommon(y,age,gender,type) else atac.glmtop.strain=edgeRfitstrain(y,age,gender,type)
    

        sel=which(atac.glmtop.strain[,"FDR"]<p.cutoff)
        write.table(cbind(annotation[ sel,"Entrez.ID"],NA,atac.glmtop.strain[sel,"logFC"]),file=paste(topgene,paste(utissue[i],int.or.slope),".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

        
    twoway.barplot.arg1=c(twoway.barplot.arg1,rep(paste(utissue[i]),2))
    twoway.barplot.arg2=c(twoway.barplot.arg2,c( sum(atac.glmtop.strain[,"FDR"]<p.cutoff & atac.glmtop.strain[,"logFC"]>0), -sum(atac.glmtop.strain[,"FDR"]<p.cutoff & atac.glmtop.strain[,"logFC"]<0)))

    twoway.barplot.arg3=c(twoway.barplot.arg3,c("+","-"))

    
#### draw heatmap of age-increasing/decreasing peaks or genes  #################
        if(length(sel)>5)
          {
    atac.glmtop=atac.glmtop.strain
    p.mat.strain=cbind(p.mat.strain,atac.glmtop[,"PValue"])
    fc.mat.strain=cbind(fc.mat.strain,atac.glmtop[,"logFC"])
    heatmapmat=y0forheatmap[atac.glmtop[,"FDR"]<p.cutoff,tissue0==utissue[i]]

    
    annot.row=annot=atac.glmtop[atac.glmtop[,"FDR"]<p.cutoff,"logFC"]

    annot[annot.row>0]=positivepeak
    annot[annot.row<0]=negativepeak
 
 
    
    colnames(heatmapmat)=paste( type,"_", samplemouseid0[tissue0==utissue[i]] ,"_",age,sep="")
   

    heatmapmat=heatmapmat[,order(age)]
    type=type[order(age)]
    heatmapmat=heatmapmat[,order(type)]


    if(nrow(heatmapmat)>40000)
      {
        tmpsample=sample(1:nrow(heatmapmat),40000)
        heatmapmat=heatmapmat[tmpsample,]
        annot=annot[tmpsample]
      }

    annot=as.data.frame(annot)
    rownames(annot)=rownames(heatmapmat)

        if(commonpattern)
          {
            colnames(annot)="common"
          }else{        
    if(typeinteraction)
      {
        colnames(annot)="slope"
      }else{
        colnames(annot)="intercept"
      }
  }
    if(min(heatmapmat)==0) heatmapmat=heatmapmat+1
   
    pheatmap(log(heatmapmat),scale="row",cluster_cols = FALSE,main=paste(utissue[i], sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]>0),positivepeak, sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]<0),negativepeak),annotation_row=annot,show_rownames=F,color=colorRampPalette(c("blue","white","red"))(100))  
 #   
colnames(fc.mat.strain)[ncol(fc.mat.strain)]=colnames(p.mat.strain)[ncol(p.mat.strain)]=paste(utissue[i],int.or.slope)
  }
      }
  }
save(p.mat.strain,fc.mat.strain,file=paste("pmat_fcmat_",int.or.slope,".Rdata",sep=""))


if(tid<5)
  {
#### heatmap of p-values of peaks/genes across tissue 
global.heatmap(p.mat.strain,fc.mat.strain)

### peaks/genes that are commonly increasing/decreasing across tissues
common.peaks(p.mat.strain,fc.mat.strain,FALSE,topgene,annotation)


all.increasing=nrow(read.delim(paste(topgene,"_all_increasing",".txt",sep="")))
all.decreasing=nrow(read.delim(paste(topgene,"_all_decreasing",".txt",sep="")))


twoway.barplot.arg1=c(twoway.barplot.arg1,rep("common",2))
twoway.barplot.arg2=c(twoway.barplot.arg2,c(all.increasing,-all.decreasing))
twoway.barplot.arg3=c(twoway.barplot.arg3,c("+","-"))

}

ylimmax=max(abs(twoway.barplot.arg2))
ylimmax=40000;YLIM=c(-ylimmax,ylimmax)

q1=twoway.barplot(twoway.barplot.arg1,twoway.barplot.arg2,twoway.barplot.arg3,(-10):10*1000,(-10):10*1000,"Tissue","no. differential peaks/genes",paste("B6 vs NZO",int.or.slope),YLIM)


multiplot(q1,NA,NA,NA,cols=2)



##

tissue.gender.type <- colnames(p.mat.strain)

if(tid<5)
  {
### save differential peaks/genes of each tissue as a txt file ####

p.mat=p.mat.strain
fc.mat=fc.mat.strain
#diff.peaks(p.mat,fc.mat,topgene)
  
### See if the age-related pattern is common across tissues


    for(jj in 1:2)
      {
        if(jj==1) tt=((p.mat<0.001 & fc.mat>0)) else tt=((p.mat<0.001 & fc.mat<0))

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
        print(mean(abs(fisher.p-fisher.p0),na.rm=T)) ## to check if my calculation is correct
        fisher.p=signif(fisher.p,2)
        fisher.p[upper.tri(fisher.p,diag=T)]="*"
        fisher.stat[upper.tri(fisher.stat,diag=T)]= 0

        if(jj==1)
          {
            write.csv(fisher.p,file=paste("fisher_pvalue_increasing_strain.csv",sep=""),quote=F)
            pheatmap(fisher.stat,scale="none",cluster_cols = FALSE,cluster_rows=FALSE,main=paste("overlap of peaks",positivepeak,"(odds ratio)"))

          } else
        {
          write.csv(fisher.p,file=paste("fisher_pvalue_decreasing_strain.csv",sep=""),quote=F)
          pheatmap(fisher.stat,scale="none",cluster_cols = FALSE,cluster_rows=FALSE,main=paste("overlap of peaks",negativepeak,"(odds ratio)"))

        }
      }
}
dev.off()  
#### Do pathway analysis using immune genes ####

library("biomaRt")

load("../../ATACseq/data/biomaRt_human_mouse.Rdata")

load("../../ATACseq/data/mousehumangene_annotation.Rdata")

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
    temptissue=tissue.gender.type[k]

    if(paste(topgene,temptissue,".txt",sep="") %in% dir())  scanfile=scan(paste(topgene,temptissue,".txt",sep=""))

    if(length(scanfile)>2)
      {
   
        diff.gene=as.matrix(read.table(paste(topgene,temptissue,".txt",sep=""),header=F))### differential gene
        if(directionsel==2)
          {
            diff.gene=rbind(diff.gene[diff.gene[,3]>0,])
            directionselsymbol="positive"
          }
        
        if(directionsel==3)
          {
            diff.gene=rbind(diff.gene[diff.gene[,3]<0,])
            directionselsymbol="negative"
          }
        
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
        write.table(human.diff.gene,file=paste(topgene,temptissue,directionselsymbol,"human.txt",sep=""),quote=F,row.names=F,col.names=F)

        if(pathwaytype==2)
          {
            setwd("~/homer")
            
            write.table(human.diff.gene,file=paste("temp",tid,".txt",sep=""),quote=F,row.names=F,col.names=F)

            system(paste("findGO.pl temp",tid,".txt human temp",tid,sep=""))
            wiki=as.matrix(read.delim(paste("temp",tid,"/wikipathways.txt",sep="")))


            wiki=rbind(wiki[as.numeric(wiki[,3])< 10^(-4) ,c(2,3)])
       # wiki=wiki[match(unique(wiki[,1]),wiki[,1]),]
        if(nrow(wiki)>0) enrichpath.wiki[[N]][[k]]=wiki

        wiki=as.matrix(read.delim(paste("temp",tid,"/kegg.txt",sep="")))
        wiki=rbind(wiki[as.numeric(wiki[,3])< 10^(-4) ,c(2,3)])
        if(nrow(wiki)>0)
          {
        wiki=wiki[match(unique(wiki[,1]),wiki[,1]),]
        enrichpath.kegg[[N]][[k]]=wiki
      }
        setwd(previous.dir)
      }
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
                
#                if(pathp[i]<0.01)
#                  {
#                    print(paste("all pathways including unnamed with p<0.01 for",temptissue))
#                    print(path.annotation[ match(pathid[i],path.annotation[,1]), 2])
#                    print( cbind(genesV2[temp,2],gene.and.CA[match(genesV2[  temp  ,1],gene.and.CA[,1]),3]))
#                  }
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
all.path.res[[3]]=enrichpath.wiki
all.path.res[[4]]=enrichpath.kegg

 save(all.path.res,file="enrichpathwayStrain.Rdata")


### draw plots of pathway analysis results ###

if(commonpattern)
  {

    pdf(paste("pathway",int.or.slope,"pattern.pdf",sep="_"))

  }else{
if(typeinteraction) pdf(file="pathwayplot_B6_vs_NZO_slope.pdf") else  pdf(file="pathwayplot_B6_vs_NZO_intercept.pdf")
  }

source("../../ATACseq/src/plot.r")


dev.off()
