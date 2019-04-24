#### This code does pathway enrichment analysis of the gene set "human.diff.gene" #####
balloonplot.pathway.cutoff=0.1
if(pathwaytype==2)
  {
    setwd("~/homer")
    source("homer_pathway_enrichment.r")                        
    
    for(pathwayid in 1:length(enrichpath.wiki))
      if(nrow(wiki[[pathwayid]])>0) enrichpath.wiki[[pathwayid]][[N]][[k]]=wiki[[pathwayid]]
    
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
if(celltype.annotation2)
  {
    load("../../ATACseq/data/geneset.info.RData")
    all.gene=  geneset.genes.scrnaseq_pbmc_specific
    path.annotation= geneset.names.scrnaseq_pbmc_specific
  }
if(celltype.annotation3)
  {
    load("../../ATACseq/data/geneset.info.RData")
    all.gene=  geneset.genes.scrnaseq_pbmc_simple_specific
    path.annotation= geneset.names.scrnaseq_pbmc_simple_specific
  }
if(celltype.annotation4)
  {
    load("../../ATACseq/data/geneset.info.RData")
    all.gene=  geneset.genes.dice_cd4
    path.annotation= geneset.names.dice_cd4
  }
if(celltype.annotation5)
  {
    load("../../ATACseq/data/geneset.info.RData")
    all.gene=  geneset.genes.dice_major
    path.annotation= geneset.names.dice_major
  }
if(celltype.annotation6)
  {
    load("../../ATACseq/data/geneset.info.RData")
    all.gene=  geneset.genes.dice_pooled
    path.annotation= geneset.names.dice_pooled
  }

tmp=paste(path.annotation[,2]," (",add.number(all.gene,path.annotation),")",sep="")
path.annotation[,2]=tmp  


module.enrich=module.enrichment.test(human.diff.gene,all.gene,gene.universe)
pathid=module.enrich[[1]];
pathp=module.enrich[[2]]
pathno=module.enrich[[3]]
genesinmodule=module.enrich[[4]]

allpath=rbind(allpath,cbind(path.annotation[ match(pathid,path.annotation[,1]), 2],pathp,pathno,length(human.diff.gene)))
                                        #allpath=rbind(allpath,c("NA",NA,sum(human.diff.gene %in% all.gene)))                    

pick=rep(NA,nrow(allpath))
for(i in 1:length(pick))
  pick[i]=pick_null_pathway(allpath[i,1])
allpath=allpath[!pick,]

enrichpath[[N]][[k]]=rbind(allpath[p.adjust(as.numeric(allpath[,2]),"fdr")<balloonplot.pathway.cutoff & allpath[,1]!="Unknown",1:2])


nopath[[N]][[k]]=rbind(allpath[ allpath[,1]!="Unknown",c(1,3,4)])
