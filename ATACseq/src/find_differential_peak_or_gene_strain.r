
p.cutoff=0.01

p.mat.strain=NULL ##bed[,"p.value"]
fc.mat.strain= NULL  ##-bed[,"Fold"]

if(tid==3) typeinteraction=TRUE else typeinteraction=FALSE

#### remove samples of age 26.75 months ######

selectsample= (AGE != 26.75)
y0=Y[,selectsample]
age0=AGE[selectsample]
gender0=GENDER[selectsample]
tissue0=TISSUE[selectsample]
type0=TYPE[selectsample]
samplemouseid0=SAMPLEMOUSEID[selectsample]

## Do differential analysis between NZO and B6 in each tissue type ######

setwd("../results")
utissue=c( "spleen", "BM", "memory","naive", "PBL")         

if(typeinteraction) pdf("diffpeakstrain_B6_vs_NZO_slope.pdf") else  pdf("diffpeakstrain_B6_vs_NZO_intercept.pdf")

topgene="topgene"

for(i in 1:length(utissue))
  {
    y=y0[,tissue0==utissue[i]]
    age=age0[tissue0==utissue[i]]
    gender=gender0[tissue0==utissue[i]]
    type=type0[tissue0==utissue[i]]
    if(length(unique(age[type=="B6"]))>1 & length(unique(age[type=="NZO"]))>1)
      {
        ### Do differential analysis between NZO and B6
    atac.glmtop.strain=edgeRfitstrain(y,age,gender,type)
   
    print("peaks significantly different between strains")
    print( sum(atac.glmtop.strain[,"FDR"]<p.cutoff))

#### draw heatmap of age-increasing/decreasing peaks or genes  #################
    
    atac.glmtop=atac.glmtop.strain
    p.mat.strain=cbind(p.mat.strain,atac.glmtop[,"PValue"])
    fc.mat.strain=cbind(fc.mat.strain,atac.glmtop[,"logFC"])
    heatmapmat=y0[atac.glmtop[,"FDR"]<p.cutoff,tissue0==utissue[i]]$counts

    
    annot.row=annot=atac.glmtop[atac.glmtop[,"FDR"]<p.cutoff,"logFC"]
    annot[annot.row>0]="higher in NZO"
    annot[annot.row<0]="lower in NZO"

    
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

    if(typeinteraction)
      {
        colnames(annot)="slope"
      }else{
        colnames(annot)="intercept"
      }

    if(min(heatmapmat)==0) heatmapmat=heatmapmat+1
   
    pheatmap(log(heatmapmat),scale="row",cluster_cols = FALSE,main=paste(utissue[i], sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]>0),"higher", sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]<0),"lower"),annotation_row=annot,show_rownames=F)
    
 #   
colnames(fc.mat.strain)[ncol(fc.mat.strain)]=colnames(p.mat.strain)[ncol(p.mat.strain)]=paste(utissue[i],"strain")
  }
  }

#### heatmap of p-values of peaks/genes across tissue 
global.heatmap(p.mat.strain,fc.mat.strain)

dev.off()

### peaks/genes that are commonly increasing/decreasing across tissues
common.peaks(p.mat.strain,fc.mat.strain,FALSE,topgene,annotation)



tissue.gender.type <- colnames(p.mat.strain)


### save differential peaks/genes of each tissue as a txt file ####

p.mat=p.mat.strain
fc.mat=fc.mat.strain
diff.peaks(p.mat,fc.mat,topgene)

  

### See if the age-related pattern is common across tissues


    for(jj in 1:2)
      {
        if(jj==1) tt=((p.mat<0.001 & fc.mat>0)) else tt=((p.mat<0.001 & fc.mat<0))

        fisher.p=fisher.p0=matrix(NA,nr=ncol(p.mat),nc=ncol(p.mat))
        rownames(fisher.p)=colnames(fisher.p)=colnames(p.mat)

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
                }
            }
        fisher.p-fisher.p0 ## to check if my calculation is correct
        fisher.p=signif(fisher.p,2)
        fisher.p[upper.tri(fisher.p,diag=T)]="*"

        if(jj==1) write.csv(fisher.p,file=paste("fisher_pvalue_increasing_strain.csv",sep=""),quote=F) else write.csv(fisher.p,file=paste("fisher_pvalue_decreasing_strain.csv",sep=""),quote=F)

      }
  
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
 save(all.path.res,file="enrichpathwayStrain.Rdata")


### draw plots of pathway analysis results ###
if(typeinteraction) pdf(file="pathwayplot_B6_vs_NZO_slope.pdf") else  pdf(file="pathwayplot_B6_vs_NZO_intercept.pdf")
source("../../ATACseq/src/plot.r")


dev.off()
