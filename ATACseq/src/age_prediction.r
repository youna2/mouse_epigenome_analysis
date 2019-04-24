library(glmnet)
source("../src/library_function.r")
remove.CD8 <- function(x)
  {
    for(i in 1:length(x))
      {
        tmp=strsplit(x[i],"CD8.")[[1]]
        x[i]=tmp[length(tmp)]
      }
    return(x)    

  }
predict.res <- function(pred_elastic,AGE,STRAIN,TISSUE,GENDER,title)
  {
    df=data.frame(x=AGE,y=pred_elastic,STRAIN=STRAIN,TISSUE=remove.CD8(TISSUE))
    temp=cor.test(df$x,df$y)
    mad=abs(df$x-df$y)
    MIN=min(c(df$x,df$y),na.rm=T)
    MAX=max(c(df$x,df$y),na.rm=T)
    
    y1 <- ggplot(df,aes(x,y))+ geom_point( aes(color=STRAIN)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = (MIN+MAX)/2, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)#+ stat_smooth(method="lm",col="grey", se=FALSE)


    y2 <- ggplot(df,aes(x,y))+ geom_point( aes(color=TISSUE)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = (MIN+MAX)/2, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)+theme(plot.title = element_text(size = 12))#+ scale_colour_manual(  values=c("PBL"="red","spleen"= "blue"))#+ stat_smooth(method="lm",col="grey", se=FALSE)

    y3 <- ggplot(df,aes(x,y))+ geom_point( aes(color=GENDER)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = (MIN+MAX)/2, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)#+ stat_smooth(method="lm",col="grey", se=FALSE)

    multiplot(y1,NA,y2,NA,y3,NA,cols=3)

#    print(c("lasso",cor.test(pred_lasso,AGE)$estimate,mean(abs(pred_lasso-AGE))))
    print(c("elastic",cor.test(pred_elastic,AGE)$estimate,mean(abs(pred_elastic-AGE))))
#    print(c("ridge",cor.test(pred_ridge,AGE)$estimate,mean(abs(pred_ridge-AGE))))

    return(list(mad,y2))
}

CVfit <- function(bed,AGE,TISSUE,SAMPLEMOUSEID)
{
  bed0=bed
  AGE0=AGE
  SAMPLEMOUSEID0=SAMPLEMOUSEID
  
utissue=unique(TISSUE)

for(k in 1:length(utissue))
  {
    bed=bed0[,TISSUE==utissue[k]]
    AGE=AGE0[TISSUE==utissue[k]]
    SAMPLEMOUSEID=SAMPLEMOUSEID0[TISSUE==utissue[k]]
    uSAMPLEMOUSEID=unique(SAMPLEMOUSEID)
    pred_elastic=pred_ridge=pred_lasso=rep(NA,length(SAMPLEMOUSEID))
    coef.est=NULL
    
    for(j in 1:length(uSAMPLEMOUSEID))
      {
        sel=which(SAMPLEMOUSEID==uSAMPLEMOUSEID[j])
        x_training=t(bed)[-sel,]
        y_training=AGE[-sel]

        x_testing=rbind(t(bed)[sel,])
        y_testing=AGE[sel]


        fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
        fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
        fit_ridge <- cv.glmnet(x=x_training, y=y_training,  alpha=0)

#        coef.est=rbind(coef.est,coef(fit_elastic)[-1,])
        ### this option performs worse
#        pred_lasso0[j] <- predict(fit_lasso, s=fit_lasso$lambda.1se, newx=x_testing)  
#        pred_elastic0[j] <- predict(fit_elastic, s=fit_lasso$lambda.1se, newx=x_testing)
#        pred_ridge0[j] <- predict(fit_ridge, s=fit_lasso$lambda.1se, newx=x_testing)


        pred_lasso[sel] <- predict(fit_lasso, s="lambda.min", newx=x_testing)
print("sel");print(sel);print(pred_elastic)
        pred_elastic[sel] <- predict(fit_elastic, s="lambda.min", newx=x_testing)
        pred_ridge[sel] <- predict(fit_ridge, s="lambda.min", newx=x_testing)

      }
    save(pred_lasso,pred_ridge,pred_elastic,file=paste("CV_",utissue[k],".Rdata",sep="") )
  }

}


# If you want to regenerate the bed file, for ATAC-seq you need to rerun ../../ATACseq/src/runall_analysis.r  until remove_problemsample.r with BTID = 6, if not, you can load the bed file
# For RNAseq, to regenerate the bed file, you need to rerun ../../RNAseq/src/runall_analysis.r until remove_problemsample.r, If not, you can load bed file###
if(ATAC) load("../../ATACseq/data/ATACseqData_for_age_prediction.Rdata")
if(RNA) load("../../RNAseq/data/RNAseqData_for_age_prediction.Rdata")

rownames(bed)=annotation.orig[,"Entrez.ID"]
bed=bed[!is.na(annotation.orig[,"Entrez.ID"]),]
annotation.orig=cbind(annotation.orig[!is.na(annotation.orig[,"Entrez.ID"]),])


### For flow data, to regenerate the bed file, you need to rerun ../../flow/src/runall_analysis.r until before TRANSFORM=F, if not, you can load bed file ###

if(flow) load("../../flow/data/FlowData_for_age_prediction.Rdata")

### This is for human data #######
### For mouse, set humandata =FALSE, for human, humandata=TRUE,
### For RNA, set RNA=TRUE. For ATAC, set ATAC=TRUE
if(humandata)
  {
    if(RNA)
      {
        bed=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_rna_glm.results.txt"))
        rownames(bed)=bed[,"GeneName"]    
        bed0=bed[,1:43]    
      }

    if(ATAC)
      {
        bed=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/pbmc_whitelisted_filtered_glm_filtered.txt")) ### This is the human ATAC file after filitering se
        rownames(bed)=bed[,"GeneName"]    
        bed0=bed[,11:54]
      }
    
    bed=matrix(as.numeric(bed0),nr=nrow(bed0))
    colnames(bed)=colnames(bed0)
    rownames(bed)=rownames(bed0)
    mean(as.numeric(bed0)==bed)
    
    sample=as.matrix(read.delim("../../ATACseq/data/aging_analysis_human_Eladio/firstcohort/sample.info.comprehensive.pbmc.txt"))
    AGE=as.numeric(sample[,"age"][match(colnames(bed),sample[,1])])
    GENDER=sample[,"sex"][match(colnames(bed),sample[,1])]
    TISSUE=rep("PBMC",length(AGE))
    STRAIN=rep("human",length(AGE))
    SAMPLEMOUSEID=colnames(bed)
  }



### This is to see how using only the genes in each immune module affect the prediction accuracy
source("../../ATACseq/src/generate_gene_list_per_immune_module.r")

### prediction using all genes and all samples

CVfit(bed,AGE,rep("alldata",length(AGE)),SAMPLEMOUSEID)

### prediction using all genes in each tissue type
CVfit(bed,AGE,TISSUE,SAMPLEMOUSEID)

## #CVfit(bed,AGE,STRAIN,SAMPLEMOUSEID)
## #CVfit(bed,AGE,GENDER,SAMPLEMOUSEID)

##This is to compare FC as well as mean expression of inflammation genes between B6 and NZO##

source("../../ATACseq/src/immune_module_genes_plot.r")


### This is to draw absolute difference of accuracy across tissues or strains
pdf(file="prediction_using_all_genes.pdf",width=10)

load(paste("CV_alldata.Rdata"))

alltissue=predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"all tissue")[[1]]
alltissueplot=predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"all tissue")[[2]]
#predict.res(pred_ridge,AGE,STRAIN,TISSUE,GENDER,"all tissue")
#predict.res(pred_lasso,AGE,STRAIN,TISSUE,GENDER,"all tissue")
names(alltissue)=rep(paste("all tissue (",length(alltissue),")",sep=""),length(alltissue))

utissue=sort(unique(TISSUE))
bytissue=NULL
for(k in 1:length(utissue))
  {
    load(paste("CV_",utissue[k],".Rdata",sep="") )
    bytissuek=predict.res(pred_elastic,AGE[TISSUE==utissue[k]],STRAIN[TISSUE==utissue[k]],TISSUE[TISSUE==utissue[k]],GENDER[TISSUE==utissue[k]],utissue[k])[[1]]
    names(bytissuek)=rep(paste(utissue[k]," (",length(bytissuek),")",sep=""),length(bytissuek))
    bytissue=c(bytissue,bytissuek)
 #   predict.res(pred_ridge,AGE[TISSUE==utissue[k]],STRAIN[TISSUE==utissue[k]],TISSUE[TISSUE==utissue[k]],GENDER[TISSUE==utissue[k]],utissue[k])
 #   predict.res(pred_lasso,AGE[TISSUE==utissue[k]],STRAIN[TISSUE==utissue[k]],TISSUE[TISSUE==utissue[k]],GENDER[TISSUE==utissue[k]],utissue[k])
    
  }

## ustrain=sort(unique(STRAIN))
## bystrain=rep(NA,length(ustrain));names(bystrain)=ustrain
## for(k in 1:length(ustrain))
##   {
##     load(paste("CV_",ustrain[k],".Rdata",sep="") )
##     bystrain[k]=predict.res(pred_elastic,AGE[STRAIN==ustrain[k]],STRAIN[STRAIN==ustrain[k]],TISSUE[STRAIN==ustrain[k]],GENDER[STRAIN==ustrain[k]],ustrain[k])
##   }

## ugender=sort(unique(GENDER))
## bygender=rep(NA,length(ugender));names(bygender)=ugender
## for(k in 1:length(ugender))
##   {
##     load(paste("CV_",ugender[k],".Rdata",sep="") )
##     bygender[k]=predict.res(pred_elastic,AGE[GENDER==ugender[k]],STRAIN[GENDER==ugender[k]],TISSUE[GENDER==ugender[k]],GENDER[GENDER==ugender[k]],ugender[k])
##   }

df=data.frame(x=c(names(alltissue),names(bytissue)),y=c(alltissue,bytissue))
save(df,alltissueplot,file="prediction_accuracy.Rdata")

fontsize=11

p <- ggplot(df, aes(x, y)) + geom_boxplot()+labs(x="tissue",y="| predicted-true ages | in months",title="Accuracy of prediction from ATACseq")+ theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ scale_fill_manual(values=c( "gold"))


 p


### This is to get the predictor 

bed=bed-rowMeans(bed)
sdvec=apply(bed,1,sd)
sdvec[sdvec==0]=1
x_training=t(bed/sdvec)
rownames(x_training)=1:nrow(x_training)
y_training=AGE

fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
fit_ridge <- cv.glmnet(x=x_training, y=y_training,  alpha=0)

#plot(fit_lasso)
#coef(fit_lasso)
#colnames(x_training)[coef(fit_lasso)[-1,]!=0]   
title=""


coef.est=coef(fit_elastic)[-1,]

COR.Cutoff=0.8
COR=cor(x_training,x_training[,coef.est>0])
extended.positive.predictor.id=colnames(x_training)[which(COR>COR.Cutoff,arr.ind=T)[,"row"]]

COR=cor(x_training,x_training[,coef.est<0])
extended.negative.predictor.id=colnames(x_training)[which(COR>COR.Cutoff,arr.ind=T)[,"row"]]

emat=x_training[,coef.est!=0]
predictor=coef.est[coef.est!=0]
print(mean(names(predictor)==colnames(emat)))
names(predictor)=convert.mouse.entrez.to.human.symbol(colnames(emat))

extended.positive.predictor.id=convert.mouse.entrez.to.human.symbol(extended.positive.predictor.id)
extended.negative.predictor.id=convert.mouse.entrez.to.human.symbol(extended.negative.predictor.id)

mean(names(predictor)[predictor>0] %in% extended.positive.predictor.id)
mean(names(predictor)[predictor<0] %in% extended.negative.predictor.id)

save(predictor,extended.positive.predictor.id,extended.negative.predictor.id,file="predictor.Rdata")

colnames(emat)=paste(1:ncol(emat),predictor)
coef.est=coef.est[coef.est!=0]

ord <- order(TISSUE,AGE,STRAIN)
emat <- emat[ord,]
anno.row <- data.frame(tissue=TISSUE[ord],age=AGE[ord], strain=STRAIN[ord])
rownames(anno.row) <- rownames(emat)

ord <- order(coef.est)
emat <- emat[,ord]
anno.col <- data.frame(coef.positive=factor(coef.est[ord]>0))
rownames(anno.col) <- colnames(emat)
#ann_colors=list(coef=c("blue","white","red"))

pheatmap(emat, annotation_row=anno.row, annotation_col=anno.col,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=F,main=title,scale="column")
dev.off()
