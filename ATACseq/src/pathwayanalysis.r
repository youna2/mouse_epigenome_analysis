### This code does the pathway enrichment analysis of human.diff.gene for immune modules, kegg and wiki pathway
library("biomaRt")

load("../../ATACseq/data/biomaRt_human_mouse.Rdata")

load("../../ATACseq/data/mousehumangene_annotation.Rdata")
gene.universe=toupper(unique(humangene[,1])) ### This defines universe of genes used for enrichment test. If doing enrichment test for mouse, the gene universe should be homologous genes (which is a saved object in the above Rdata file), and this command should not be used. 


setwd("../results")
previous.dir=getwd()

fdr.cutoff=0.1
all.path.res=all.path.res2=vector("list",4)
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

    enrichpath=nopath=enrichpath.wiki=enrichpath.kegg=vector("list",3)
    for(i in 1:3) enrichpath[[i]]=nopath[[i]]=enrichpath.wiki[[i]]=enrichpath.kegg[[i]]=vector("list",length(tissue.gender.type))
    
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
    all.path.res2[[pathwaytype]]=nopath
  }
all.path.res[[3]]=enrichpath.wiki
all.path.res[[4]]=enrichpath.kegg
