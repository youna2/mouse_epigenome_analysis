### This code does the pathway enrichment analysis of human.diff.gene for immune modules, kegg and wiki pathway

source("../../ATACseq/src/library_function.r")

library("biomaRt")

load("../../ATACseq/data/biomaRt_human_mouse.Rdata")

load("../../ATACseq/data/mousehumangene_annotation.Rdata")
#gene.universe=toupper(unique(humangene[,1])) ### This defines universe of genes used for enrichment test. If doing enrichment test for mouse, the gene universe should be homologous genes (which is a saved object in the above Rdata file), and this command should not be used. 


setwd("../results")
previous.dir=getwd()


tid=1
k=1
N=1

#human.diff.gene=names(predictor)[!is.na(names(predictor))]



fdr.cutoff=0.1





all.path.res=all.path.res2=vector("list",25)

for(pathwaytype in 1:7)
  {
     immunemodule=celltype.annotation=celltype.annotation2=celltype.annotation3=celltype.annotation4=celltype.annotation5=celltype.annotation6=FALSE
    if(pathwaytype==1)
      celltype.annotation=TRUE
      
    if(pathwaytype==2)
      immunemodule=TRUE

    if(pathwaytype==3)
      celltype.annotation2=TRUE
  
    if(pathwaytype==4)
      celltype.annotation3=TRUE
  
    if(pathwaytype==5)
      celltype.annotation4=TRUE
  
    if(pathwaytype==6)
      celltype.annotation5=TRUE

    if(pathwaytype==7)
      celltype.annotation6=TRUE

    if(pathwaytype==2)
      {
        enrichpath.wiki=vector("list",18)
        for(count in 1:length(enrichpath.wiki))
          {
            enrichpath.wiki[[count]]=vector("list",3)
            for(i in 1:3)
              enrichpath.wiki[[count]][[i]]=vector("list",1)
          }
      }
    
    enrichpath=nopath=vector("list",3)
    for(i in 1:3) enrichpath[[i]]=nopath[[i]]=vector("list",1)


    source("../../ATACseq/src/pathway_enrichment_test.r")
                
    all.path.res[[pathwaytype]]=enrichpath
    all.path.res2[[pathwaytype]]=nopath
  }

for(count in 1:length(enrichpath.wiki))
all.path.res[[count+7]]=enrichpath.wiki[[count]]



