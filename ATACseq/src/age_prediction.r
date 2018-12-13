library(glmnet)
predict.res <- function(pred_elastic,AGE,STRAIN,TISSUE,GENDER,title)
  {
    df=data.frame(x=AGE,y=pred_elastic,STRAIN=STRAIN,TISSUE=TISSUE)
    temp=cor.test(df$x,df$y)
    mad=abs(df$x-df$y)
    MIN=min(c(df$x,df$y),na.rm=T)
    MAX=max(c(df$x,df$y),na.rm=T)
    
    y1 <- ggplot(df,aes(x,y))+ geom_point( aes(color=STRAIN)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = 5+MIN, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)#+ stat_smooth(method="lm",col="grey", se=FALSE)


    y2 <- ggplot(df,aes(x,y))+ geom_point( aes(color=TISSUE)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle("(A) Predicted age from flow data")+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = 5+MIN, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)+theme(plot.title = element_text(size = 12))+ scale_colour_manual(  values=c("PBL"="red","spleen"= "blue"))#+ stat_smooth(method="lm",col="grey", se=FALSE)

    y3 <- ggplot(df,aes(x,y))+ geom_point( aes(color=GENDER)) +ylab(paste("predicted age"))+xlab(paste("true age"))+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2)), x = 5+MIN, y = MAX, size = 6, colour = "red")+ geom_abline(intercept=0,slope=1)+xlim(MIN,MAX)+ylim(MIN,MAX)#+ stat_smooth(method="lm",col="grey", se=FALSE)

    multiplot(y1,y2,y3,NA,cols=2)

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


#        fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
        fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
#        fit_ridge <- cv.glmnet(x=x_training, y=y_training,  alpha=0)

        coef.est=rbind(coef.est,coef(fit_elastic)[-1,])
        ### this option performs worse
#        pred_lasso0[j] <- predict(fit_lasso, s=fit_lasso$lambda.1se, newx=x_testing)  
#        pred_elastic0[j] <- predict(fit_elastic, s=fit_lasso$lambda.1se, newx=x_testing)
#        pred_ridge0[j] <- predict(fit_ridge, s=fit_lasso$lambda.1se, newx=x_testing)


#        pred_lasso[j] <- predict(fit_lasso, s="lambda.min", newx=x_testing)
print("sel");print(sel);print(pred_elastic)
        pred_elastic[sel] <- predict(fit_elastic, s="lambda.min", newx=x_testing)
#        pred_ridge[j] <- predict(fit_ridge, s="lambda.min", newx=x_testing)

      }
    save(coef.est,pred_elastic,file=paste("CV_",utissue[k],".Rdata",sep="") )
  }

}
CVfit(bed,AGE,rep("alldata",length(AGE)),SAMPLEMOUSEID)
CVfit(bed,AGE,TISSUE,SAMPLEMOUSEID)
#CVfit(bed,AGE,STRAIN,SAMPLEMOUSEID)
#CVfit(bed,AGE,GENDER,SAMPLEMOUSEID)

rownames(bed)=annotation[,"Entrez.ID"]

#pdf(file="temp.pdf",width=9)
load(paste("CV_alldata.Rdata"))

alltissue=predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"all tissue")[[1]]
alltissueplot=predict.res(pred_elastic,AGE,STRAIN,TISSUE,GENDER,"all tissue")[[2]]
names(alltissue)=rep(paste("all tissue (",length(alltissue),")",sep=""),length(alltissue))

utissue=sort(unique(TISSUE))
bytissue=NULL
for(k in 1:length(utissue))
  {
    load(paste("CV_",utissue[k],".Rdata",sep="") )
    bytissuek=predict.res(pred_elastic,AGE[TISSUE==utissue[k]],STRAIN[TISSUE==utissue[k]],TISSUE[TISSUE==utissue[k]],GENDER[TISSUE==utissue[k]],utissue[k])[[1]]
    names(bytissuek)=rep(paste(utissue[k]," (",length(bytissuek),")",sep=""),length(bytissuek))
    bytissue=c(bytissue,bytissuek)
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

## df$x=factor(df$x,levels=df$x)
## p<-ggplot(data=df, aes(x=x, y=y,fill=class)) + geom_bar(stat="identity")+ylab("median absoute error in months")
## p

sdvec=apply(bed,1,sd)
sdvec[sdvec==0]=1
x_training=t(bed/sdvec)
rownames(x_training)=1:nrow(x_training)
y_training=AGE
#fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)


#plot(fit_lasso)
#coef(fit_lasso)
#colnames(x_training)[coef(fit_lasso)[-1,]!=0]   
title=""


coef.est=coef(fit_elastic)[-1,]
emat=x_training[,coef.est!=0]

predictor=coef.est[coef.est!=0]
names(predictor)=colnames(emat)
names(predictor)=convert.mouse.entrez.to.human.symbol(colnames(emat))
save(predictor,file="predictor.Rdata")

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
#dev.off()
