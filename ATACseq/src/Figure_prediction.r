library(ggplot2)
source("../../ATACseq/src/library_function.r")
before.paren <- function(x)
  {
    for(i in 1:length(x))
      x[i]=strsplit(x[i]," [(]")[[1]][1]
    return(x)
  }
fontsize=11
load("../../ATACseq/src/prediction_accuracy.Rdata")
ATAC=paste("ATAC (",sort(table(df$x),decreasing=T)[1],")",sep="")
df$z=rep(ATAC,nrow(df))
dfATAC=df

load("../../RNAseq/src/prediction_accuracy.Rdata")
RNA=paste("RNA (",sort(table(df$x),decreasing=T)[1],")",sep="")
df$z=rep(RNA,nrow(df))
dfRNA=df

load("../../flow/src/prediction_accuracy.Rdata")
flow=paste("flow (",sort(table(df$x),decreasing=T)[1],")",sep="")
df$z=rep(flow,nrow(df))
dfFLOW=df
alltissueplotflow=alltissueplot

df0=rbind(dfATAC,dfRNA,dfFLOW)

df0$z=factor(df0$z,levels=c(flow,ATAC,RNA),ordered=TRUE)

df=df0[before.paren(as.character(df0$x))=="all tissue",]


p1 <- ggplot(df, aes(z, y,fill=z)) + geom_boxplot()+labs(x="Data",y="| predicted-true ages | in months",title="(C) Age prediction using all tissues")+ylim(0,10)+theme(plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ scale_fill_manual(values=c("indianred2", "gold", "seagreen2"))



df=df0[df0$z==flow,]

p2 <- ggplot(df, aes(x, y,fill=z)) + geom_boxplot()+labs(x="tissue",y="| predicted-true ages | in months",title="(B) Accuracy of prediction from Flow")+ylim(0,10)+ theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ scale_fill_manual(values=c("indianred2"))



df=df0[df0$z==ATAC,]

p3 <- ggplot(df, aes(x, y,fill=z)) + geom_boxplot()+labs(x="tissue",y="| predicted-true ages | in months",title="(D) Accuracy of prediction from ATACseq")+ylim(0,10)+ theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ scale_fill_manual(values=c( "gold"))

df=df0[df0$z==RNA,]

p4 <- ggplot(df, aes(x, y,fill=z)) + geom_boxplot()+labs(x="tissue",y="| predicted-true ages | in months",title="(E) Accuracy of prediction from RNAseq")+ylim(0,10)+ theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = fontsize,hjust=0),legend.position="none")+ scale_fill_manual(values=c("seagreen2"))




load("../../flow/src/predictor.Rdata")
names(predictor)[names(predictor)=="X..Viable"]="Viable cell percentage"
names(predictor)[names(predictor)=="CD4.Lymphocytes"]="CD4"
df=data.frame(x=names(predictor),y=predictor)
p<-ggplot(data=df, aes(x=reorder(x,y), y=y,fill=y>0)) + geom_bar(stat="identity")+labs(x="predictors",y="coefficient",title="(C) Selected predictors from Flow")+coord_flip()+theme(plot.title = element_text(size = fontsize,hjust=1),legend.position="none")+ scale_fill_manual(values=c("darkorchid1","violetred"))

pdf(file="Figure_prediction.pdf",width=10)
multiplot(alltissueplotflow,p3,p2,p4,p,NA,cols=3)
dev.off()
