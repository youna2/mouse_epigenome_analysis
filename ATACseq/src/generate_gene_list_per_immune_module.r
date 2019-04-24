RN=100 
version="2008"#"2015"
all.gene=as.matrix(read.table(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_genes.txt",sep=""),header=T))
path.annotation=as.matrix(read.csv(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_annotations.csv",sep="")))
inflammation1=all.gene[all.gene[,1]=="M3_2",2]
save(inflammation1,file="../results/inflammation1_gene_list.Rdata")
Tcell=all.gene[all.gene[,1]=="M2_8",2]
save(Tcell,file="../results/tcell_gene_list.Rdata")


madvector=madp=NULL
pdf(file="prediction_using_genes_in_each_module.pdf")

CVfit(bed,AGE,rep("alldata",length(AGE)),SAMPLEMOUSEID)
load(paste("CV_alldata.Rdata"))
madvector=rbind(madvector,cbind("all genes",(predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"all genes")[[1]])))

set.seed(10)

#predict.res(pred_ridge,AGE,STRAIN,TISSUE,GENDER,"all tissue")
#predict.res(pred_lasso,AGE,STRAIN,TISSUE,GENDER,"all tissue")



for(i in 1:nrow(path.annotation))
  {
    if(strsplit(path.annotation[i,2],"")[[1]][1]!="U")
    {
      module=all.gene[all.gene[,1]==path.annotation[i,1],2]

      if(humandata) sel=toupper(rownames(bed)) %in% toupper(module) else sel=toupper(convert.mouse.entrez.to.human.symbol(rownames(bed))) %in% toupper(module)
      CVfit(bed[sel,],AGE,rep("alldata",length(AGE)),SAMPLEMOUSEID)
      load(paste("CV_alldata.Rdata"))

      tmpmodule=predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,strsplit(path.annotation[i,2],"/")[[1]][1])[[1]]
      madvector=rbind(madvector,cbind(strsplit(path.annotation[i,2],"/")[[1]][1],tmpmodule))

      madrandom=rep(NA,RN)
      for(ii in 1:RN)
        {
          if(humandata) sel=toupper(rownames(bed)) %in% toupper(module) else sel=toupper(convert.mouse.entrez.to.human.symbol(rownames(bed))) %in% toupper(module)
          CVfit(bed[sample((1:nrow(bed))[!sel],sum(sel)),],AGE,rep("alldata",length(AGE)),SAMPLEMOUSEID)
          load(paste("CV_alldata.Rdata"))
          madrandom[ii]=mean(predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"Random")[[1]])
        }

      madp=rbind(madp,c(strsplit(path.annotation[i,2],"/")[[1]][1],mean(madrandom<=mean(tmpmodule))))
      
#predict.res(pred_ridge,AGE,STRAIN,TISSUE,GENDER,"all tissue")
#predict.res(pred_lasso,AGE,STRAIN,TISSUE,GENDER,"all tissue")
    }
  }
fontsize=11
df=data.frame(x=madvector[,1],y=as.numeric(madvector[,2]),z= (madvector[,1]=="all genes" | madvector[,1]=="Inflammation I"))

ggplot(df, aes(x, y,fill=z)) + geom_boxplot()+labs(x="Immune module",y="| predicted-true ages | in months",title="Age prediction using all tissues")+theme(plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=c("indianred2", "gold", "seagreen2"))

dev.off()
save(madp,file="prediction_pvalue_immune_module.Rdata")
