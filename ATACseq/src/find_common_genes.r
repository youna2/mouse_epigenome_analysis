### For RNA

gene.and.annotation <- function(x,all.gene)
  {
    x=as.vector(x)
    tmp=matrix(0,nr=length(x),nc=length(unique(all.gene[,1])))
    colnames(tmp)=unique(all.gene[,1])
    for(i in 1:length(x))
      {
        tmp[i,all.gene[toupper(all.gene[,2])==toupper(x[i]),1]]=1
      }
    return(cbind(x,tmp))
  }
convert.to.path <- function(x,path.annotation)
  {
    for(i in 1:nrow(path.annotation))
      x[x==path.annotation[i,1]]=path.annotation[i,2]
    return(x)
  }

intersect.matrix <- function(x,y)
  {
    tmp=intersect(x[,1],y[,1])
    return(x[x[,1]%in% tmp,])
  }

load("../../ATACseq/data/pbmc_specific_genes.annotations_April.2018.EJM_10x.RData")
all.gene=  geneset.genes.scrnaseq_pbmc_specific
path.annotation= geneset.names.scrnaseq_pbmc_specific[,c(2,1)]
all.gene[,1]=convert.to.path(all.gene[,1],path.annotation)


TISSUE=c("PBL","spleen","CD8.naive","CD8.memory","PBL.PBMC.common")
no=c(3,2,5,4)
STRAIN=c("NZO","B6")
DIRECTION=c("positive","negative")

a=b=vector("list",2)
RNAlist=ATAClist=RNA.ATAC.list=vector("list",length(DIRECTION))
names(RNAlist)=names(ATAClist)=names(RNA.ATAC.list)=DIRECTION


for(i in 1:length(DIRECTION))
  {
    direction=DIRECTION[i]

    RNAlist[[i]]=ATAClist[[i]]=RNA.ATAC.list[[i]]=vector("list",length(TISSUE))
    names(RNAlist[[i]])=names(ATAClist[[i]])=names(RNA.ATAC.list[[i]])=TISSUE

    for(k in 1:(length(TISSUE)-1))
      {
        tissue=TISSUE[k]
        RNAlist[[i]][[k]]=ATAClist[[i]][[k]]=RNA.ATAC.list[[i]][[k]]=vector("list",length(STRAIN))
        names(RNAlist[[i]][[k]])=names(ATAClist[[i]][[k]])=names(RNA.ATAC.list[[i]][[k]])=STRAIN

        for(j in 1:length(STRAIN))
          {
            strain=STRAIN[j]
            RNAlist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../RNAseq/results/within ",strain,tissue," within ",strain,direction,"human.txt",sep=""))),all.gene)
            ATAClist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../ATACseq/results",no[k],"/within ",strain,tissue," within ",strain,direction,"human.txt",sep=""))),all.gene)
            
            RNA.ATAC.list[[i]][[k]][[j]]=intersect.matrix(RNAlist[[i]][[k]][[j]],ATAClist[[i]][[k]][[j]])

            print(c(nrow(RNAlist[[i]][[k]][[j]]),nrow(ATAClist[[i]][[k]][[j]]),nrow(RNA.ATAC.list[[i]][[k]][[j]])))
          }
      }



    k=length(TISSUE)
    tissue=TISSUE[k]
    if(i==1) dir="increasing" else dir="decreasing"
    RNAlist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../RNAseq/results/common_",direction,"_genes.txt",sep=""))),all.gene)
    ATAClist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../ATACseq/results/common_",direction,"_genes.txt",sep=""))),all.gene)
            
    RNA.ATAC.list[[i]][[k]][[j]]=intersect.matrix(RNAlist[[i]][[k]][[j]],ATAClist[[i]][[k]][[j]])


    
  }

gene=ATAClist$positive$PBL$B6
gene[gene[,"T cells"]==1,1]


## tmp1=tmp2=tmp3=NULL
## for(k in 1:length(TISSUE))
## {
##   tmp1=c(tmp1,RNAlist[[i]][[k]])
##   tmp2=c(tmp2,ATAClist[[i]][[k]])
##   tmp3=c(tmp3,RNA.ATAC.list[[i]][[k]])
## }
