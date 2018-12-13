library(ggplot2)
library(MASS)
library(viridis)
theme_set(theme_bw(base_size = 16))
#### generating plot of FC in B6 vs FC in NZO #####

source("../../ATACseq/src/library_function.r")


immunemodule=TRUE
celltype.annotation=FALSE


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



for(RNA in c(TRUE,FALSE))
  {
    if(RNA) load("../../RNAseq/results/RNAseqData2.Rdata") else load("../../ATACseq/results/ATACseqData2.Rdata")

    dim(annotation)
    
    pathid=unique(all.gene[,1])
    pathway=matrix(NA,nr=nrow(annotation),nc=2+length(pathid))
    pathway[,1]=annotation[,"Entrez.ID"]
    pathway[,2]=convert.mouse.entrez.to.human.symbol(annotation[,"Entrez.ID"])



    for(i in 1:length(pathid))
      {
        path.gene=toupper(all.gene[all.gene[,1]==pathid[i],2])### all genes in pathway i
        pathway[,i+2]= pathway[,2] %in% path.gene
      }


    mean(path.annotation[ match(pathid,path.annotation[,1]), 1]==pathid)

    colnames(pathway)=c("mouse_entrezid","human_gene_symbol",path.annotation[ match(pathid,path.annotation[,1]), 2])



    if(RNA) load("../../RNAseq/results/pmat_fcmat_within B6.Rdata") else load("../../ATACseq/results/pmat_fcmat_within B6.Rdata")
    fc.B6=fc.mat.strain
    p.B6=p.mat.strain
    if(RNA) load("../../RNAseq/results/pmat_fcmat_within NZO.Rdata") else load("../../ATACseq/results/pmat_fcmat_within NZO.Rdata")
    fc.NZO=fc.mat.strain
    p.NZO=p.mat.strain
    if(RNA)
      {
        pdf(file="FC_B6_vs_NZO_RNA.pdf",width=9)
        datatype="RNA";max=0.6;textsize=6
      } else
    {
      pdf(file="FC_B6_vs_NZO_ATAC.pdf",width=9)
      datatype="ATAC";max=0.2;textsize=6
    }

    for(N in 1:3)
      {
        if(N==1)
          {
            B6significant=FALSE
            NZOsignificant=FALSE
          }

        if(N==2)
          {
            B6significant=TRUE
            NZOsignificant=FALSE
          }

        if(N==3)
          {
            B6significant=FALSE
            NZOsignificant=TRUE
          }

        total.FC.graph=vector("list",ncol(fc.B6))    
        allpath=path.annotation[c(1:3,5,7,9:12,14,16,20:22),2]
        corvec=rep(NA,ncol(fc.B6))
        for(i in 1:ncol(fc.B6))
          {
            count=j=0
            corvec[i]=cor.test(fc.B6[,i],fc.NZO[,i])$estimate

            df=data.frame(B6=fc.B6[,i],NZO=fc.NZO[,i])
            df$density <- get_density(df$B6,df$NZO)
            
            total.FC.graph[[i]] <- ggplot(df)+ geom_point( aes(B6, NZO,color=density),size=1/10) + scale_color_viridis() +xlab(paste("log FC in",strsplit(colnames(fc.B6)[i]," ")[[1]][3]))+ylab(paste("log FC in",strsplit(colnames(fc.NZO)[i]," ")[[1]][3]))+ggtitle(strsplit(colnames(fc.B6)[i]," ")[[1]][1])+geom_abline(color="grey",slope=1,intercept=0)+xlim(-max,max)+ylim(-max,max)+annotate("text", label = paste("r=",signif(corvec[i],2)), x = 0, y = max, size = textsize, colour = "red")+ theme(legend.position="none")




            
            while(j <= length(allpath))
              {
                p=vector("list",4)
                
                for(k in 1:4)
                  {
                    j=k+count*4

                    if(j <= length(allpath))
                      {
                        selpath=allpath[j]
                        df=data.frame(B6=fc.B6[,i],NZO=fc.NZO[,i],annot=pathway[,selpath],significanceB6=p.adjust(p.B6[,i],"fdr")<0.05,significanceNZO=p.adjust(p.NZO[,i],"fdr")<0.05)
                        df=df[order(df$annot),]
                        
                        if(B6significant) df$B6[!df$significanceB6]=NA
                        if(NZOsignificant) df$NZO[!df$significanceNZO]=NA
                        
                        p[[k]] <- ggplot(df, aes(B6, NZO))+ geom_point(aes(colour=annot,size=annot))+ scale_color_manual(values=c( "FALSE"="blue","TRUE"="red"))+ scale_size_manual( values=c("FALSE"=1/10,"TRUE"=1/2))+xlab(colnames(fc.B6)[i])+ylab(colnames(fc.NZO)[i])+ggtitle(selpath)+geom_abline(color="grey",slope=1,intercept=0)
                      }
                  }

                count=count+1
                multiplot(p[[1]],p[[2]],p[[3]],p[[4]],cols=2)  
              }
          }
      }
    dev.off()



    if(RNA)
      {
        total.FC.graph.RNA=total.FC.graph
        save(corvec,file="../../ATACseq/results/cor_B6_NZO_RNA.Rdata")

      } else{
        total.FC.graph.ATC=total.FC.graph
        save(corvec,file="../../ATACseq/results/cor_B6_NZO_ATAC.Rdata")
      }
  }

draw.heatmap=T
if(draw.heatmap)
  {
    library(pheatmap)
    load("../../ATACseq/results/cor_B6_NZO_RNA.Rdata")
    corvecRNA=corvec
    load("../../ATACseq/results/cor_B6_NZO_ATAC.Rdata")
    corvecATAC=corvec
    cormat=rbind(corvecRNA,corvecATAC)
    colnames(cormat)=as.matrix(as.data.frame(strsplit(colnames(fc.B6)," ")))[1,]
    rownames(cormat)=c("RNA","ATAC")
    pheatmap(cormat,cluster_cols = FALSE,cluster_rows = FALSE,scale="none",main="Correlation of fold change between B6 and NZO",color = colorRampPalette(c("blue","white","red"))(100))
    
    pdf(file="total_FC_graph.pdf",width=20)
    multiplot(total.FC.graph.ATC[[1]],total.FC.graph.RNA[[1]],total.FC.graph.ATC[[2]],total.FC.graph.RNA[[2]],total.FC.graph.ATC[[3]],total.FC.graph.RNA[[3]],total.FC.graph.ATC[[4]],total.FC.graph.RNA[[4]],total.FC.graph.ATC[[5]],total.FC.graph.RNA[[5]],cols=5)
    dev.off()
  }
