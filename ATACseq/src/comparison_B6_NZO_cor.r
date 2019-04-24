
   version="2008"#"2015"
all.gene=as.matrix(read.table(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_genes.txt",sep=""),header=T))
path.annotation=as.matrix(read.csv(paste("../../ATACseq/data/immunemodule/VP",version,"_Modules_annotations.csv",sep="")))

   strain.cor.module=NULL
   tmp=rep(NA,ncol(B6FC))
   for(k in 1:ncol(B6FC)) tmp[k]=cor(B6FC[,k],NZOFC[,k])
   strain.cor.module=rbind(strain.cor.module,c("allgene",tmp))

for(i in 1:nrow(path.annotation))
  {
   if(strsplit(path.annotation[i,2],"")[[1]][1]!="U")
     {
       tmp=rep(NA,ncol(B6FC))
      module=all.gene[all.gene[,1]==path.annotation[i,1],2]
      sel=toupper(humangenename) %in% toupper(module)
      for(k in 1:ncol(B6FC)) tmp[k]=cor(B6FC[sel,k],NZOFC[sel,k])
       strain.cor.module=rbind(strain.cor.module,c(path.annotation[i,2],tmp))
     }
 }


cbind(strain.cor.module[,1],rowMeans(matrix(as.numeric(strain.cor.module[,-1]),nrow=nrow(strain.cor.module))))

   df=data.frame(module=rep(strain.cor.module[,1],4),tissue=rep(convert.col.name(colnames(B6FC)),each=nrow(strain.cor.module)),cor=as.numeric(strain.cor.module[,-1]))

   p1=ggplot(data=df, aes(x=module, y=cor)) +  geom_bar(stat="identity", position=position_dodge())+  theme_minimal()+ggtitle("B6 vs. NZO")+coord_flip()+labs(x="")+facet_wrap(~tissue,nrow=1)
   
   
