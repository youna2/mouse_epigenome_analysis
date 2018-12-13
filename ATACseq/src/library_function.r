library(ggplot2)
library(MASS)
library(viridis)
library(gtable)
library(gridExtra)
plot.correlationplot <- function(x,y,title,xlab,ylab)
  { 
    df=data.frame(ATAC=x,RNA=y)
    df$density <- get_density(df$ATAC,df$RNA)
    if(RNA) margin=10^(-3) else margin=0
    
    temp=cor.test(x,y)
    y <- ggplot(df,aes(ATAC, RNA))+ geom_point( aes(ATAC, RNA,color=density),size=1/10) + scale_color_viridis() +xlab(xlab)+ylab(ylab)+ggtitle(title)+annotate("text", label = paste("r=",signif(temp$estimate,2),sep=""), x = 0, y = quantile(y,1-margin), size = textsize, colour = "red")+ theme(legend.position="none") + stat_smooth(method="lm",col="grey", se=FALSE)+xlim(quantile(x,margin),quantile(x,1-margin))+ylim(quantile(y,margin),quantile(y,1-margin))

    return(y)
  }


plot.balloonplot.mouse <- function(balloonplot)
  {

    for(i in 1:3)
      {
        multiballoonplot=vector("list",4)
    
        for(j in 1:3)
          {
            pathwaymat=balloonplot[[i]]

            if(i==3)
              {
                pathwaymat=pathwaymat[as.numeric(pathwaymat[,2])<10^(-5),]
                pathwaymat=pathwaymat[pathwaymat[,4]=="CD8.memory" | pathwaymat[,4]=="CD8.naive" ,]
              }
        
            if(i!=3) pathwaymat=pathwaymat[pathwaymat[,4]!="CD8.memory" &pathwaymat[,4]!="CD8.naive" ,]
            
            if(j ==1) sel=pathwaymat[,5]=="within NZO" | pathwaymat[,5]=="within B6"

            if(j ==2) sel=pathwaymat[,5]=="B6 vs NZO common" | pathwaymat[,5]=="B6 vs NZO opposite"
            if(j ==3) sel=pathwaymat[,5]=="only NZO significant" | pathwaymat[,5]=="only B6 significant"

            if(j ==4) sel=pathwaymat[,5]=="intercept" | pathwaymat[,5]=="slope"        

            pathwaymat=pathwaymat[sel,]
    
            df=data.frame(pathid=pathwaymat[,1],logp=-log10(as.numeric(pathwaymat[,2])+10^(-30)),direction=(pathwaymat[,3]),tissue=pathwaymat[,4],class=pathwaymat[,5])

            multiballoonplot[[j]]=(ggballoonplot(df, x = "tissue", y = "pathid", size = "logp", fill = "direction" ,facet.by="class",ggtheme = theme_bw()))+scale_fill_manual(values = c("down"="dodgerblue","up"="firebrick1"),guide=guide_legend(title = NULL))# +facet_wrap("class",scales="free_x") #+
#  scale_fill_viridis_c(option = "C")
          }

#        if(i<3)
#          {
            multiplot(multiballoonplot[[1]],multiballoonplot[[2]],multiballoonplot[[3]],multiballoonplot[[4]],col=1)
      ##     }else{
      ##       multiplot(multiballoonplot[[1]],multiballoonplot[[2]],col=1);
      ##       print(multiballoonplot[[4]]);

      ## }
  }

}

plot.balloonplot.human <- function(balloonplot)
  {
    multiballoonplot=vector("list",3)
    for(i in 1:3)
      {    
        pathwaymat=balloonplot[[i]]
        
        pathwaymat=pathwaymat[pathwaymat[,4]=="PBL" | pathwaymat[,5]=="human",]

        if(i==3)
          {
            pathwaymat=pathwaymat[as.numeric(pathwaymat[,2])<10^(-3),]
          }
                    
        sel=pathwaymat[,5]=="within NZO" | pathwaymat[,5]=="within B6"| pathwaymat[,5]=="human"


        pathwaymat=pathwaymat[sel,]
    
        df=data.frame(pathid=pathwaymat[,1],logp=-log10(as.numeric(pathwaymat[,2])+10^(-30)),direction=(pathwaymat[,3]),tissue=pathwaymat[,4],class=pathwaymat[,5])

        multiballoonplot[[i]]=(ggballoonplot(df, x = "tissue", y = "pathid", size = "logp", fill = "direction",facet.by="class" ,ggtheme = theme_bw()))+scale_fill_manual(values = c("down"="dodgerblue","up"="firebrick1"),guide=guide_legend(title = NULL))# +facet_wrap("class",scales="free_x") #+
#  scale_fill_viridis_c(option = "C")
      }

    p1=multiballoonplot[[1]]
    p2=multiballoonplot[[2]]
    p3=multiballoonplot[[3]]
    grid.arrange(p1,p2,p3,nrow=3)

}

module.enrichment.test <- function(human.diff.gene,all.gene,gene.universe)
  {
print(c("geneuniverse",mean(toupper(human.diff.gene) %in% gene.universe)))
    human.diff.gene=intersect(toupper(human.diff.gene),gene.universe)
                    pathid=unique(all.gene[,1])
                    pathp=rep(NA,length(pathid))
                    for(i in 1:length(pathid))
                      {
                        path.gene=toupper(all.gene[all.gene[,1]==pathid[i],2])### all genes in pathway i

                        total=length(unique(gene.universe))  ### all human ortholog genes 
                        white=length(unique(intersect(path.gene,gene.universe))) ### all genes in pathway i in universe
                        black=total-white
                        pick=length(human.diff.gene)

                        intersection=intersect(human.diff.gene,path.gene)
                        
                        lintersection= length(intersection )

                        whitepick=lintersection-1

                        pathp[i]= 1-phyper(whitepick,white,black,pick)

#                        temp=match(intersection,genesV2[,2])
                        
                      }

                    return(list(pathid,pathp))
  }


## convert.mouse.entrez.to.human.symbol <- function(x)
##   {
##     library("biomaRt")
##     load("../../ATACseq/data/biomaRt_human_mouse.Rdata")
##     load("../../ATACseq/data/mousehumangene_annotation.Rdata")

##     genesV2 = getLDS(attributes = c("entrezgene"), filters = "entrezgene", values =as.numeric(x) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

##     genesV2[,2]=toupper(genesV2[, 2])
##     gene.and.CA=genesV2[match(as.numeric(x),genesV2[,1]),]

##     print(mean(gene.and.CA[,1]==as.numeric(x),na.rm=T))
##     return(gene.and.CA[,2])
##   }
###

convert.mouse.entrez.to.human.symbol <- function(x)
  {
    genesV2=read.csv("../../ATACseq/data/Human.Mouse.Rat.Ortholog.Gene.List.csv",header=T)
    genesV2=genesV2[,c("Mouse_Entrez","Human_Symbol")]
    
    genesV2[,2]=toupper(genesV2[, 2])
    gene.and.CA=genesV2[match(as.numeric(x),genesV2[,1]),]

    print(mean(gene.and.CA[,1]==as.numeric(x),na.rm=T))
    return(gene.and.CA[,2])
  }


###
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


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
#                scale_y_continuous(labels = c(4,2,0,2,4)*10000,breaks = c(-4,-2,0,2,4)*10000) +
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
df=NULL
    for(i in 1:ncol(p.mat))
      {
        temptissue=colnames(p.mat)[i]
        print(temptissue)
        print(sum(p.mat[,i]<0.05))
        sel=which(p.mat[,i]<0.05)
        write.table(cbind(annotation[ sel,"Entrez.ID"],p.mat[sel,i],fc.mat[sel,i]),file=paste(topgene,temptissue,".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")


        
        sel=which(p.mat[,i]<0.05 & fc.mat[,i]>0)
        if(length(sel)>0) df=rbind(df,   cbind( paste(colnames(p.mat)[i],"+"), removing.paren(annotation[sel,"Annotation"])))
        sel=which(p.mat[,i]<0.05 & fc.mat[,i]<0)
        if(length(sel)>0) df=rbind(df,   cbind(paste(colnames(p.mat)[i],"-"),removing.paren(annotation[sel,"Annotation"])))

      }
temp=table(df[,1],df[,2])
temp=temp/rowSums(temp)
df=data.frame(tissue=rep(rownames(temp),ncol(temp)),percent=c(temp),annotation=rep(colnames(temp),each=nrow(temp)))

p <- ggplot(data=df, aes(x=tissue, y=percent, fill=annotation)) +
  geom_bar(stat="identity")+   scale_fill_brewer(palette = "Set3") +   ylab("Percent") +    ggtitle("Homer annotation of peaks")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

print( p)
        
      
  }
removing.paren <- function(x)
{
  y=x
  for(i in 1:length(x)) y[i]=strsplit(x[i],"[(]")[[1]][1]
  return(y)
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

edgeRfitcommon <- function(y,age,gender,type)
  {

    if(tid==5)
      {
    sel= (type=="B6")
    resB6=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    sel= (type=="NZO")
    resNZO=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])

    res=matrix(1,nr=nrow(resB6),nc=ncol(resB6))
    colnames(res)=colnames(resB6)
    res[,"logFC"]=0
    res[resB6[,"logFC"]<0 & resNZO[,"logFC"]<0,"logFC"]= -1
    res[resB6[,"logFC"]>0 & resNZO[,"logFC"]>0,"logFC"]= 1

    res[,"FDR"]=1
    res[resB6[,"FDR"]< p.cutoff & resNZO[,"FDR"]<p.cutoff & res[,"logFC"]!=0 ,"FDR"]= p.cutoff/2
  }

 if(tid==6)
      {
    sel= (type=="B6")
    resB6=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    sel= (type=="NZO")
    resNZO=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])

    res=matrix(1,nr=nrow(resB6),nc=ncol(resB6))
    colnames(res)=colnames(resB6)
    res[,"logFC"]=0
    res[resB6[,"logFC"]>0 & resNZO[,"logFC"]<0,"logFC"]= -1
    res[resB6[,"logFC"]<0 & resNZO[,"logFC"]>0,"logFC"]= 1

    res[,"FDR"]=1
    res[resB6[,"FDR"]< p.cutoff & resNZO[,"FDR"]<p.cutoff & res[,"logFC"]!=0 ,"FDR"]= p.cutoff/2
  }

     if(tid==7)
      {
    sel= (type=="B6")
    resB6=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    sel= (type=="NZO")
    resNZO=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])

    res=matrix(1,nr=nrow(resB6),nc=ncol(resB6))
    colnames(res)=colnames(resB6)
    res[,"logFC"]=0
    res[ resNZO[,"logFC"]<0,"logFC"]= -1
    res[ resNZO[,"logFC"]>0,"logFC"]= 1

    res[,"FDR"]=1
    res[resB6[,"FDR"]>= p.cutoff & resNZO[,"FDR"]<p.cutoff & res[,"logFC"]!=0 ,"FDR"]= p.cutoff/2
  }

     if(tid==8)
      {
    sel= (type=="B6")
    resB6=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    sel= (type=="NZO")
    resNZO=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])

    res=matrix(1,nr=nrow(resB6),nc=ncol(resB6))
    colnames(res)=colnames(resB6)
    res[,"logFC"]=0
    res[ resB6[,"logFC"]<0,"logFC"]= -1
    res[ resB6[,"logFC"]>0,"logFC"]= 1

    res[,"FDR"]=1
    res[resNZO[,"FDR"]>= p.cutoff & resB6[,"FDR"]<p.cutoff & res[,"logFC"]!=0 ,"FDR"]= p.cutoff/2
  }




###
    if(tid==9)
      {
        sel= (type=="B6")
        resB6=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
        res=resB6
      }
    if(tid==10)
      {
        sel= (type=="NZO")
        resNZO=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
        res=resNZO
      }
    if(tid==11)
      res=edgeRfit(y,age,gender,type)
    
    ## sel= (type=="B6" & gender =="M")
    ## resB6M=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    ## sel= (type=="NZO" & gender=="M")
    ## resNZOM=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])

    ## sel= (type=="B6" & gender =="F")
    ## resB6F=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
    ## sel= (type=="NZO" & gender=="F")
    ## resNZOF=edgeRfit(y[,sel],age[sel],gender[sel],type[sel])
 

    ## res=matrix(1,nr=nrow(resB6M),nc=ncol(resB6M))
    ## colnames(res)=colnames(resB6M)
    ## res[,"logFC"]=0
    ## res[resB6M[,"logFC"]<0 & resNZOM[,"logFC"]<0 & resB6F[,"logFC"]<0 & resNZOF[,"logFC"]<0,"logFC"]= -1
    ## res[resB6M[,"logFC"]>0 & resNZOM[,"logFC"]>0 & resB6F[,"logFC"]>0 & resNZOF[,"logFC"]>0,"logFC"]= 1

    ## res[,"FDR"]=1
    ## res[resB6M[,"FDR"]< p.cutoff & resNZOM[,"FDR"]<p.cutoff & resB6F[,"FDR"]< p.cutoff & resNZOF[,"FDR"]<p.cutoff & res[,"logFC"]!=0 ,"FDR"]= p.cutoff/2


    return(res)    
    
  }


    

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
#    adj.p.mat=p.mat
#    for(i in 1:ncol(p.mat)) adj.p.mat[,i]=p.adjust(p.mat[,i],"fdr")
    x=x[rowSums(p.mat<p.cutoff)>0,]

    if(nrow(x)>40000) x=x[sample(1:nrow(x),40000),]
    if(max(x)> abs(min(x))) x=rbind(x,-max(abs(x)))  else  x=rbind(x,max(abs(x)))#for color balance
    
    if(tid==1 | tid==2 | tid==5)
      {
        pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log q* sign of change" )
      }
    else
      {
        if(typeinteraction)
          {
            pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log q* sign of change of slope in NZO" )        
          }
        else
          {
            pheatmap(x,scale="none",cluster_cols = FALSE, color = colorRampPalette(c("red4", "white", "blue4"))(100),show_rownames = F,main="log q* sign of change of intercept in NZO" )        
          }
      }
    

    
  }

