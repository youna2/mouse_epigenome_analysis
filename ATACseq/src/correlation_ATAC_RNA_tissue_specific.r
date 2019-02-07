library(ggplot2)
library(MASS)
library(viridis)
theme_set(theme_bw(base_size = 16))

### generating plot of FC in RNA vs FC in ATAC ####

source("../../ATACseq/src/library_function.r")
RNA=FALSE

textsize=6
for(selB6 in c(TRUE,FALSE))
  {
    corvec=rep(NA,5)
    total.FC.graph=vector("list",5)    

    for(BTID in 2:5)
      {

        load("../../RNAseq/results/RNAseqData2.Rdata")


        if(selB6) load("../../RNAseq/results/pmat_fcmat_within B6.Rdata") else load("../../RNAseq/results/pmat_fcmat_within NZO.Rdata")
        p.mat=p.mat.strain
        fc.mat=fc.mat.strain
        ## removeNOFC=F
        ## if(removeNOFC)
        ##   {
        ##     sel=(abs(fc.mat[,i])>0.0001)
        ##     annotation=cbind(annotation[sel,])
        ##     colnames(annotation)="Entrez.ID"
        ##     fc.mat=fc.mat[sel,]
        ##     p.mat=p.mat[sel,]
        ##   }


        EntrezRNA=annotation[,"Entrez.ID"]
        annotationRNA=annotation
        p.mat.RNA=p.mat
        fc.mat.RNA=fc.mat

        load(paste("../../ATACseq/results",BTID,"/ATACseqData2.Rdata",sep=""))

        if(selB6) load(paste("../../ATACseq/results",BTID,"/pmat_fcmat_within B6.Rdata",sep="")) else load(paste("../../ATACseq/results",BTID,"/pmat_fcmat_within NZO.Rdata",sep=""))
        p.mat=p.mat.strain
        fc.mat=fc.mat.strain


        onlypromoter=F ### this doesn't make much difference
        if(onlypromoter)
          {
            sel=(removing.paren(annotation[,"Annotation"])=="promoter-TSS ")
            annotation=annotation[sel,]
            fc.mat=cbind(fc.mat[sel,])
            p.mat=cbind(p.mat[sel,])
          }

        EntrezATAC=annotation[,"Entrez.ID"]
        annotationATAC=annotation
        p.mat.ATAC=p.mat
        fc.mat.ATAC=fc.mat

        p.mat.RNA=cbind(p.mat.RNA[,colnames(p.mat.ATAC)])
        fc.mat.RNA=cbind(fc.mat.RNA[,colnames(fc.mat.ATAC)])
        
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

                                        #geom_abline(color="grey",slope=1,intercept=0)


        ATAC.FC=cbind( -sign(fc.mat.ATAC[matching.index[,1],])*log(p.mat.ATAC[matching.index[,1],]))
        RNA.FC=cbind( -sign(fc.mat.RNA[matching.index[,2],])*log(p.mat.RNA[matching.index[,2],]))

        mouseentrez=common


        for(i in 1:ncol(fc.mat.ATAC))
          {
            x=fc.mat.ATAC[matching.index[,1],i]
            y=fc.mat.RNA[matching.index[,2],i]

            x.p=p.mat.ATAC[matching.index[,1],i]
            y.p=p.mat.RNA[matching.index[,2],i]

                                        #        x=x[x.p<0.1]
                                        #        y=y[x.p<0.1]
            
            title=strsplit(colnames(fc.mat.ATAC)[i]," ")[[1]][1]
            xlab="log FC in ATAC"
            if(i==1)
              {
                if(selB6) ylab="B6\nlog FC in RNA" else ylab="NZO\nlog FC in RNA"
              }else{
                ylab="log FC in RNA"
              }
            total.FC.graph[[BTID]] <- plot.correlationplot(x,y,title,xlab,ylab)
            
                                        #        corvec[i]=temp$estimate
          }
      }
    if(selB6)
      {
        total.FC.graph.B6=total.FC.graph
        save(ATAC.FC,RNA.FC,mouseentrez,file="../../ATACseq/results/FC_ATAC_RNA_B6.Rdata")
      }else{
        total.FC.graph.NZO=total.FC.graph
        save(ATAC.FC,RNA.FC,mouseentrez,file="../../ATACseq/results/FC_ATAC_RNA_NZO.Rdata")
      }
  }

draw.heatmap=T
if(draw.heatmap)
  {
    pdf(file="total_FC_ATAC_RNA_graph.pdf",width=20)
    multiplot(total.FC.graph.B6[[1]],total.FC.graph.NZO[[1]],total.FC.graph.B6[[2]],total.FC.graph.NZO[[2]],total.FC.graph.B6[[3]],total.FC.graph.NZO[[3]],total.FC.graph.B6[[4]],total.FC.graph.NZO[[4]],total.FC.graph.B6[[5]],total.FC.graph.NZO[[5]],cols=5)
    dev.off()
  }
