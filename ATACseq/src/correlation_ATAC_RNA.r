library(ggplot2)
library(MASS)
library(viridis)
theme_set(theme_bw(base_size = 16))

### generating plot of FC in RNA vs FC in ATAC ####

source("../../ATACseq/src/library_function.r")


textsize=6
for(selB6 in c(TRUE,FALSE))
  {
    corvec=rep(NA,5)
    total.FC.graph=vector("list",5)    
    for(i in 1:5)
      {


        load("../../RNAseq/data/RNAseqData.Rdata")


        if(selB6) load("../../RNAseq/results/pmat_fcmat_within B6.Rdata") else load("../../RNAseq/results/pmat_fcmat_within NZO.Rdata")
        p.mat=p.mat.strain
        fc.mat=fc.mat.strain
        removeNOFC=T
        if(removeNOFC)
          {
            sel=(abs(fc.mat[,i])>0.0001)
            annotation=cbind(annotation[sel,])
            colnames(annotation)="Entrez.ID"
            fc.mat=fc.mat[sel,]
            p.mat=p.mat[sel,]
          }


        EntrezRNA=annotation[,"Entrez.ID"]
        annotationRNA=annotation
        p.mat.RNA=p.mat
        fc.mat.RNA=fc.mat
        load("../../ATACseq/data/ATACseqData.Rdata")



        if(selB6) load("../../ATACseq/results/pmat_fcmat_within B6.Rdata") else load("../../ATACseq/results/pmat_fcmat_within NZO.Rdata")
        p.mat=p.mat.strain
        fc.mat=fc.mat.strain


        onlypromoter=T
        if(onlypromoter)
          {
            sel=(removing.paren(annotation[,"Annotation"])=="promoter-TSS ")
            annotation=annotation[sel,]
            fc.mat=fc.mat[sel,]
            p.mat=p.mat[sel,]
          }

        EntrezATAC=annotation[,"Entrez.ID"]
        annotationATAC=annotation
        p.mat.ATAC=p.mat
        fc.mat.ATAC=fc.mat

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
                                        #save(matching.index,file="../../ATACseq/results/matching_index.Rdata")


        df=data.frame(ATAC=fc.mat.ATAC[matching.index[,1],i],RNA=fc.mat.RNA[matching.index[,2],i])
        df$density <- get_density(df$ATAC,df$RNA)

        temp=cor.test(fc.mat.ATAC[matching.index[,1],i],fc.mat.RNA[matching.index[,2],i])
        total.FC.graph[[i]] <- ggplot(df,aes(ATAC, RNA))+ geom_point( aes(ATAC, RNA,color=density),size=1/10) + scale_color_viridis() +xlab(paste("log FC in ATAC"))+ylab(paste("log FC in RNA"))+ggtitle(strsplit(colnames(fc.mat.ATAC)[i]," ")[[1]][1])+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = 0, y = max(df$RNA,na.rm=T), size = textsize, colour = "red")+ theme(legend.position="none") + stat_smooth(method="lm",col="grey", se=FALSE)
                                        #geom_abline(color="grey",slope=1,intercept=0)
        
        corvec[i]=temp$estimate
      }
    if(selB6)
      {
        total.FC.graph.ATC=total.FC.graph
        save(corvec,file="../../ATACseq/results/cor_ATAC_RNA_B6.Rdata")
      }else{
        total.FC.graph.RNA=total.FC.graph
        save(corvec,file="../../ATACseq/results/cor_ATAC_RNA_NZO.Rdata")
      }
  }
draw.heatmap=T
if(draw.heatmap)
  {
    pdf(file="total_FC_ATAC_RNA_graph.pdf",width=20)
    multiplot(total.FC.graph.ATC[[1]],total.FC.graph.RNA[[1]],total.FC.graph.ATC[[2]],total.FC.graph.RNA[[2]],total.FC.graph.ATC[[3]],total.FC.graph.RNA[[3]],total.FC.graph.ATC[[4]],total.FC.graph.RNA[[4]],total.FC.graph.ATC[[5]],total.FC.graph.RNA[[5]],cols=5)
    dev.off()
  }
