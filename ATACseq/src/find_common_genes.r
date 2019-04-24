
genes.in.pathway <- function(x)
  {
    pathwaymat=NULL
    for(i in 1:length(x))
      {
        for(j in 1:length(x[[i]]))
          {
            if(j<5)
              {
                for(k in 1:length(x[[i]][[j]]))
                  {
                    gene=x[[i]][[j]][[k]]
                    if(!is.null(nrow(gene)))
                      for(l in select.no.U.pathway(colnames(gene)[-1]))
                        pathwaymat=rbind(pathwaymat,paste(names(x[[i]])[j],names(x[[i]][[j]])[k],names(x)[i],l,paste(gene[gene[,l]==1,1],collapse=","),sep=","))
                  }
              }else{
                gene=x[[i]][[j]]
                if(!is.null(nrow(gene)))
                  for(l in select.no.U.pathway(colnames(gene)[-1]))
                    pathwaymat=rbind(pathwaymat,paste(names(x[[i]])[j],names(x[[i]][[j]])[k],names(x)[i],l,paste(gene[gene[,l]==1,1],collapse=","),sep=","))
                
              }
          }
      }
    return(pathwaymat)
  }

gene.and.annotation <- function(x,all.gene)
  {
    x=as.vector(x)
    tmp=matrix(0,nr=length(x),nc=length(unique(all.gene[,1])))
    colnames(tmp)=unique(all.gene[,1])
    for(i in 1:length(x))
      {
        tmp[i,all.gene[toupper(all.gene[,2])==toupper(x[i]),1]]=1
      }
    return(cbind(toupper(x),tmp))
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
select.no.U.pathway <- function(x)
  {
    y=rep(T,length(x))
    for(i in 1:length(x))
      if(strsplit(x[i],"")[[1]][1]=="U") y[i]=F
    return(x[y])
  }

RNA.ATACgenelist=ATACgenelist=RNAgenelist=NULL
pathwaytitle=c("cell type","immune module","dice_pool")
for(pathwaytype in 1:3)
  {
    if(pathwaytype==1)
      {
        load("../../ATACseq/data/pbmc_specific_genes.annotations_April.2018.EJM_10x.RData")
        all.gene=  geneset.genes.scrnaseq_pbmc_specific
        path.annotation= geneset.names.scrnaseq_pbmc_specific[,c(2,1)]
        all.gene[,1]=convert.to.path(all.gene[,1],path.annotation)
      }
    if(pathwaytype==2)
      {
        version="2008"#"2015"
        all.gene=as.matrix(read.table(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_genes.txt",sep=""),header=T))
        path.annotation=as.matrix(read.csv(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_annotations.csv",sep="")))
        all.gene[,1]=convert.to.path(all.gene[,1],path.annotation)
      }
    if(pathwaytype==3)
      {
        load("../../ATACseq/data/geneset.info.RData")
        all.gene=  geneset.genes.dice_pooled
        path.annotation= geneset.names.dice_pooled
        all.gene[,1]=convert.to.path(all.gene[,1],path.annotation)
      }

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
### list of mouse genes increasing/decreasing for each tissue and strain
                strain=STRAIN[j]
                RNAlist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../RNAseq/results/within ",strain,tissue," within ",strain,direction,"human.txt",sep=""))),all.gene)
                ATAClist[[i]][[k]][[j]]=gene.and.annotation(as.matrix(read.table(paste("../../ATACseq/results",no[k],"/within ",strain,tissue," within ",strain,direction,"human.txt",sep=""))),all.gene)
                
                RNA.ATAC.list[[i]][[k]][[j]]=intersect.matrix(RNAlist[[i]][[k]][[j]],ATAClist[[i]][[k]][[j]])

                print(c(nrow(RNAlist[[i]][[k]][[j]]),nrow(ATAClist[[i]][[k]][[j]]),nrow(RNA.ATAC.list[[i]][[k]][[j]])))
              }
          }


#### list of genes commonly increasing/decreasing between human PBML and mouse PBL B6/NZO
        k=length(TISSUE)
        tissue=TISSUE[k]
        if(i==1) dir="increasing" else dir="decreasing"
        RNAlist[[i]][[k]]=gene.and.annotation(as.matrix(read.table(paste("../../RNAseq/results/common_",dir,"_genes.txt",sep=""))),all.gene)
        ATAClist[[i]][[k]]=gene.and.annotation(as.matrix(read.table(paste("../../ATACseq/results/common_",dir,"_genes.txt",sep=""))),all.gene)
        
        RNA.ATAC.list[[i]][[k]]=intersect.matrix(RNAlist[[i]][[k]],ATAClist[[i]][[k]])
        
      }
    save(RNAlist,ATAClist,RNA.ATAC.list,file="../results/list_of_common_genes.Rdata")

                                        #gene=ATAClist$negative$PBL.PBMC.common
                                        #gene=RNA.ATAC.list$negative$PBL.PBMC.common

    RNA.ATACgenelist=rbind(RNA.ATACgenelist,"","",pathwaytitle[pathwaytype],"","",genes.in.pathway(RNA.ATAC.list))
    ATACgenelist=rbind(ATACgenelist,"","",pathwaytitle[pathwaytype],"","",genes.in.pathway(ATAClist))
    RNAgenelist=rbind(RNAgenelist,"","",pathwaytitle[pathwaytype],"","",genes.in.pathway(RNAlist))
  }
write.csv(RNA.ATACgenelist,file="../results/list_of_changing_genes_RNA_ATAC.csv",quote=F,row.names=F)
write.csv(ATACgenelist,file="../results/list_of_changing_genes_ATAC.csv",quote=F,row.names=F)
write.csv(RNAgenelist,file="../results/list_of_changing_genes_RNA.csv",quote=F,row.names=F)



generate.list.genes <- function(ATAClist)
  {
    tab=rep(NA,length(ATAClist)*length(ATAClist[[1]])+2)
    count=0
    for(dir in names(ATAClist))
      {
### list of genes increasing/decreasing for each tissue, common between B6 and NZO
        for(tissue in names(ATAClist[[dir]]))
          {
            count=count+1
            names(tab)[count]=paste(dir,tissue)
            if(tissue=="PBL.PBMC.common") tab[count]=paste(sort(ATAClist[[dir]][[tissue]][,1]),collapse=",")  else  tab[count]=paste(  sort(intersect(ATAClist[[dir]][[tissue]]$B6[,1],ATAClist[[dir]][[tissue]]$NZO[,1])),collapse=",")
          }

        count=count+1
        names(tab)[count]=paste(dir,"across tissues in mouse")

        temp=NULL
        for(i in 1:4) temp=c(temp,ATAClist[[dir]][[i]]$B6[,1],ATAClist[[dir]][[i]]$NZO[,1])
        temp=table(temp)
        print(dir)
        print(table(temp))

        ## list of genes commonly increasing/decreasing across tissues and strains
        tab[count]=paste(sort(names(temp)[temp==8]),collapse=",")
      }
    tab=cbind(tab)
    return(tab)
  }


### list of genes increasing/decreasing for each tissue, common between B6 and NZO
tab <- generate.list.genes(ATAClist)
write.csv(tab,file="../results/list_of_common_genes_ATAC.csv",quote=F,row.names=T,col.names=F)


tab <- generate.list.genes(RNAlist)
write.csv(tab,file="../results/list_of_common_genes_RNA.csv",quote=F,row.names=T,col.names=F)


tab <- generate.list.genes(RNA.ATAC.list)
write.csv(tab,file="../results/list_of_common_genes_RNA_ATAC.csv",quote=F,row.names=T,col.names=F)


## tmp1=tmp2=tmp3=NULL
## for(k in 1:length(TISSUE))
## {
##   tmp1=c(tmp1,RNAlist[[i]][[k]])
##   tmp2=c(tmp2,ATAClist[[i]][[k]])
##   tmp3=c(tmp3,RNA.ATAC.list[[i]][[k]])
## }
overlap.across.tissue <- function(ATAClist)
  {
    y=vector("list",2);
    names(y)=c("increasing","decreasing")
    
    for(k in 1:2)
      {
        tab=vector("list",4)
### list of genes increasing/decreasing for each tissue, common between B6 and NZO
        for(i in 1:4)
          {
            tab[[i]]=unique(intersect(ATAClist[[k]][[i]]$B6[,1],ATAClist[[k]][[i]]$NZO[,1]))
          }
        y[[k]]=tab
      }
    return(y)
  }

tab <- overlap.across.tissue(ATAClist)
for(i in 1:2)
  {
    print(names(tab)[i])
    print(table(table(unlist(tab[[i]])))/length(table(unlist(tab[[i]]))))
  }

tab <- overlap.across.tissue(RNAlist)
for(i in 1:2)
  {
    print(names(tab)[i])
    print(table(table(unlist(tab[[i]])))/length(table(unlist(tab[[i]]))))
  }
