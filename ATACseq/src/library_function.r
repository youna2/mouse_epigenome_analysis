


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


twoway.barplot <- function(x,y,z,labels,breaks,xlab,ylab,title,YLIM)
  {
    df=data.frame(x=x,y=y,z=z)
    
    p1=ggplot(df, aes(x,y,fill=z)) +
      geom_bar(stat = "identity", position="identity",width=0.7,size=0.5,color="white") +
        geom_hline(color="black",size=0.25,yintercept = 0) +
          geom_hline(color="firebrick4",size=0.25,linetype=2,yintercept = 0) +
            scale_fill_manual(values = c("-"="dodgerblue","+"="firebrick1"),guide=guide_legend(title = NULL)) +
              coord_flip() +
#                scale_y_continuous(labels = labels,breaks = breaks) +
                  labs(x=xlab, y=ylab) +
                    theme_minimal(base_size = 8) +
                      theme(axis.text.y = element_text(size=8),aspect.ratio = 1,panel.grid.major.y = element_line(color="honeydew2"))+ggtitle(title)+ylim(YLIM)
    return(p1)
  }

twoway.barplot.strain <- function(x,y,z,labels,breaks,xlab,ylab,title,YLIM)
  {
    df=data.frame(x=x,y=y,z=z)
    
    p1=ggplot(df, aes(x,y,fill=z)) +
      geom_bar(stat = "identity", position="identity",width=0.7,size=0.5,color="white") +
        geom_hline(color="black",size=0.25,yintercept = 0) +
          geom_hline(color="firebrick4",size=0.25,linetype=2,yintercept = 0) +
            scale_fill_manual(values = c("lower in NZO"="dodgerblue","higher in NZO"="firebrick1"),guide=guide_legend(title = NULL)) +
              coord_flip() +
#                scale_y_continuous(labels = labels,breaks = breaks) +
                  labs(x=xlab, y=ylab) +
                    theme_minimal(base_size = 8) +
                      theme(axis.text.y = element_text(size=8),aspect.ratio = 1,panel.grid.major.y = element_line(color="honeydew2"))+ggtitle(title)+ylim(YLIM)
    return(p1)
  }

convert.from.ensembl.to.entrezID.mouse <- function(x)
{

annotation <- as.matrix(getBM(attributes=c("ensembl_gene_id","entrezgene"), mart=ensembl,filters="ensembl_gene_id",values=x))

mean(annotation[match(x,annotation[,1]),1]==x,na.rm=T)

y=   annotation[match(x,annotation[,1]),2]

return(y)
}


convert.from.symbol.to.ensembl.mouse <- function(x)
{

annotation <- as.matrix(getBM(attributes=c("mgi_symbol","ensembl_gene_id"), mart=ensembl,filters="mgi_symbol",values=x))

mean(annotation[match(x,annotation[,1]),1]==x,na.rm=T)

y=   annotation[match(x,annotation[,1]),2]

return(y)
}

PCA <- function(emat,color.var)
  {
emat=emat[rowSums(emat)>0,]

pca=prcomp(t(emat),center=T,scale=T)

d <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)

xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])


dat=data.frame(PC1=as.numeric(pca$x[,1]),PC2=as.numeric(pca$x[,2]),tissue=color.var[,"TISSUE"],strain=color.var[,"STRAIN"],age=color.var[,"AGE"],gender=color.var[,"GENDER"])

for(i in 3:ncol(dat))
  {

p <- ggplot(dat, aes(PC1,PC2))+ geom_point(aes(color=dat[,i]))+labs(x=xl,y=yl,color=colnames(dat)[i])

print(p)

  }

}

common.peaks <- function(p.mat,fc.mat,sexspecific,topgene,annotation)
{
  p.cut=0.05
tissuesel=1:4*2-1

  if(sexspecific)
    {
print("In Female, peaks significantly increasing across tissues")
sel <- rowMeans(p.mat[,tissuesel]<p.cut & fc.mat[,tissuesel]>0)==1### Female increasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,1),file=paste(topgene,"_F_increasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")


print("In Female, peaks significantly decreasing across tissues")
sel <- rowMeans(p.mat[,tissuesel]<p.cut & fc.mat[,tissuesel]<0)==1### Female decreasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,-1),file=paste(topgene,"_F_decreasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")


tissuesel=1:4*2

print("In Male, peaks significantly increasing across tissues")
sel <- rowMeans(p.mat[,tissuesel]<p.cut & fc.mat[,tissuesel]>0)==1### Male increasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,1),file=paste(topgene,"_M_increasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

print("In Male, peaks significantly decreasing across tissues")
sel <- rowMeans(p.mat[,tissuesel]<p.cut & fc.mat[,tissuesel]<0)==1### Male decreasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,-1),file=paste(topgene,"_M_decreasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

}
  
print("peaks significantly increasing across tissues and genders")
sel <- rowMeans(p.mat<p.cut & fc.mat>0)==1### All increasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,1),file=paste(topgene,"_all_increasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

print("peaks significantly decreasing across tissues and genders")
sel <- rowMeans(p.mat<p.cut & fc.mat<0)==1### All decreasing across tissues
print(sum(sel))
write.table(cbind(annotation[ sel,"Entrez.ID"],p.cut,-1),file=paste(topgene,"_all_decreasing",".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

}
diff.peaks <- function(p.mat,fc.mat,topgene)
  {
### differential genes for each tissue ####

    for(i in 1:ncol(p.mat))
      {
        temptissue=colnames(p.mat)[i]
        print(temptissue)
        print(sum(p.adjust(p.mat[,i],"fdr")<0.05))
        sel=which(p.adjust(p.mat[,i],"fdr")<0.05)
        write.table(cbind(annotation[ sel,"Entrez.ID"],p.mat[sel,i],fc.mat[sel,i]),file=paste(topgene,temptissue,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
      }
  }

pick_null_pathway <- function(x)
  {
    y=FALSE
    temp=strsplit(x,"_")[[1]]
    if(temp[1]=="U") y=TRUE
    temp=strsplit(x," ")[[1]]
    if(temp[1]=="Unnamed") y=TRUE
    return(y)
  }

#################
edgeRfit <- function(y,age,gender,type)
  {
    if(length(unique(gender))==1)
      {
        if(length(unique(type))==1) design <- model.matrix(~age) else design <- model.matrix(~age+type)
      } else
    {
      if(length(unique(type))==1) design <- model.matrix(~age+gender) else design <- model.matrix(~age+gender+type)
    }#### one can treat age as a continuous variable. When continuous variable, the logFC is the log fold change per increase of age by 1, for example, when comparing logFC from 1)  age=old, young vs. 2) age=0 vs 10, I confirmed that logFC of 2) is logFC of 1) divided by 10 (the sign may be opposite).  Given the link function for negative binomial regression, it is a linear trend on the scale of the model formula and linear predictor, and this corresponds to an exponential trend on the scale of the predicted counts.
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)

                                        # Edit this to change estimation method
                                        # try both methods, usually there is not a big difference
                                        # robust settings help with noise
    test.method="glm" # ql | glm

    if (test.method=="glm") {
      atac.glm <- glmFit(y, design)
      atac.glmfit <- glmLRT(atac.glm,coef=2)
    }

    atac.glmtop <- topTags(atac.glmfit,n=nrow(bed),sort.by="none", adjust.method = "fdr")$table

    return(atac.glmtop)
  }



#################
edgeRfitstrain <- function(y,age,gender,type)
  {
    design <- model.matrix(~type+age+age:type+gender)
   
    rownames(design) <- colnames(y)
    y <- estimateDisp(y, design, robust=TRUE)

    test.method="glm" # ql | glm

    if (test.method=="glm") {
      atac.glm <- glmFit(y, design)
      if(typeinteraction)  atac.glmfit <- glmLRT(atac.glm,coef="typeNZO:age") else atac.glmfit <- glmLRT(atac.glm,coef="typeNZO")
      
    }

    atac.glmtop <- topTags(atac.glmfit,n=nrow(bed),sort.by="none", adjust.method = "fdr")$table

    return(atac.glmtop)
  }


global.heatmap <- function(p.mat,fc.mat)
  {
    temp=as.vector(p.mat)
    temp2=min(temp[temp>0])
    p.mat[p.mat==0]=temp2
    x= -log(p.mat)*sign(fc.mat)
    adj.p.mat=p.mat
    for(i in 1:ncol(p.mat)) adj.p.mat[,i]=p.adjust(p.mat[,i],"fdr")
    x=x[rowSums(adj.p.mat<p.cutoff)>0,]

    if(nrow(x)>40000) x=x[sample(1:nrow(x),40000),]
    if(max(x)> abs(min(x))) x=rbind(x,-max(abs(x)))  else  x=rbind(x,max(abs(x)))#for color balance
    
    if(tid==1 | tid==2)
      {
        pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log p* sign of change" )
      }
    else
      {
        if(typeinteraction)
          {
            pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log p* sign of change of slope in NZO" )        
          }
        else
          {
            pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log p* sign of change of intercept in NZO" )        
          }
      }
    

    
  }

