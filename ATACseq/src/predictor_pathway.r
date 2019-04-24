path.vec <- function(y)
  {
    if( !is.null(y) && nrow(y)>0 )
      {
        if(ncol(y)>2) tmp=paste(y[,1],"(",signif(as.numeric(y[,2]),3),")",y[,3]) else tmp=paste(y[,1],"(",signif(as.numeric(y[,2]),3),")")
        y=paste(tmp,collapse=",")
      } else y=""
    return(y)    
  }
pathwayvector <- function(x)
  {
    y=rep(NA,2)
    names(y)=c(",pathway,cell type",",,immune module")
    for(i in 1:2)
      y[i]= path.vec(x[[i]][[1]][[1]])
    return(y)
  }

#### Age predictors from flow cytometry data ####

load("../../flow/src/predictor.Rdata")
predictor=sort(predictor[!is.na(names(predictor))])
write.csv(cbind(signif(predictor,3)),file="Flow_predictor.csv",quote=F)

#### Age predictors from ATACseq data ####

load("../../ATACseq/src/predictor.Rdata")
tmpres=NULL
predictor=sort(predictor[!is.na(names(predictor))])
write.csv(cbind(signif(predictor,3)),file="ATACseq_predictor.csv",quote=F)

    ### enriched pathways for positive predictors
human.diff.gene=names(predictor)[predictor>0 ]
ATAC.pos.predictor=paste(sort(human.diff.gene),collapse=",")
names(ATAC.pos.predictor)="ATAC,gene,positive predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,ATAC.pos.predictor,pathwayvector(all.path.res))
all.path.res
   ### enriched pathways for extended positive predictors
human.diff.gene=extended.positive.predictor.id[!is.na(extended.positive.predictor.id)]
ATAC.extended.pos.predictor=paste(sort(human.diff.gene),collapse=",")
names(ATAC.extended.pos.predictor)=",gene,extended positive predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,ATAC.extended.pos.predictor,pathwayvector(all.path.res))
all.path.res
   ### enriched pathways for negative predictors
human.diff.gene=names(predictor)[predictor<0]
ATAC.neg.predictor=paste(sort(human.diff.gene),collapse=",")
names(ATAC.neg.predictor)=",gene,negative predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,ATAC.neg.predictor,pathwayvector(all.path.res))
all.path.res
   ### enriched pathways for extended negative predictors
human.diff.gene=extended.negative.predictor.id[!is.na(extended.negative.predictor.id)]
ATAC.extended.neg.predictor=paste(sort(human.diff.gene),collapse=",")
names(ATAC.extended.neg.predictor)=",gene,extended negative predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,ATAC.extended.neg.predictor,pathwayvector(all.path.res))
all.path.res

#### Age predictor from RNAseq data ####
load("../../RNAseq/src/predictor.Rdata")
predictor=sort(predictor[!is.na(names(predictor))])
write.csv(cbind(signif(predictor,3)),file="RNAseq_predictor.csv",quote=F)
    ### enriched pathways for positive predictors
human.diff.gene=names(predictor)[predictor>0 ]
RNA.pos.predictor=paste(sort(human.diff.gene),collapse=",")
names(RNA.pos.predictor)="RNA,gene,positive predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,RNA.pos.predictor,pathwayvector(all.path.res))
all.path.res
    ### enriched pathways for extended positive predictors
human.diff.gene=extended.positive.predictor.id[!is.na(extended.positive.predictor.id)]
RNA.extended.pos.predictor=paste(sort(human.diff.gene),collapse=",")
names(RNA.extended.pos.predictor)=",gene,extended positive predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,RNA.extended.pos.predictor,pathwayvector(all.path.res))
all.path.res
    ### enriched pathways for negative predictors
human.diff.gene=names(predictor)[predictor<0]
RNA.neg.predictor=paste(sort(human.diff.gene),collapse=",")
names(RNA.neg.predictor)=",gene,negative predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,RNA.neg.predictor,pathwayvector(all.path.res))
all.path.res
    ### enriched pathways for extended negative predictors
human.diff.gene=extended.negative.predictor.id[!is.na(extended.negative.predictor.id)]
RNA.extended.neg.predictor=paste(sort(human.diff.gene),collapse=",")
names(RNA.extended.neg.predictor)=",gene,extended negative predictor"
source("../src/pathwayanalysis_general.r")
tmpres=c(tmpres,RNA.extended.neg.predictor,pathwayvector(all.path.res))
all.path.res

tab=tmpres
write.csv(cbind(tab),file="../results/age_predictor_pathway.csv",quote=F,row.names=T)





