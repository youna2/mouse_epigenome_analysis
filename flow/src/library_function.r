coef.heatmap <- function(x,title)
  {
    MAX=max(abs(-sign(x)*log(abs(x))),na.rm=T)+1;breaksList = seq(-MAX, MAX, by = 1)

    txt=matrix("",nr=nrow(x),nc=ncol(x))
    txt[,1:2][p.adjust(abs(x[,1:2]),"fdr")<0.05]="*"
    txt[,3][p.adjust(abs(x[,3]),"fdr")<0.05]="*"
    
    pheatmap(-sign(x)*log(abs(x)),cluster_cols = FALSE,cluster_rows = FALSE,scale="none",main=title,color = colorRampPalette(c("blue","white","red"))(length(breaksList)),breaks = breaksList,display_numbers=txt)
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


PCA <- function(emat,color.var)
  {
emat=emat[rowSums(emat)>0,]

pca=prcomp(t(emat),center=T,scale=T)

d <- round(pca$sdev^2/sum(pca$sdev^2)*100, digits=1)

xl <- sprintf("PC 1: %.1f %%", d[1])
yl <- sprintf("PC 2: %.1f %%", d[2])


dat=data.frame(PC1=as.numeric(pca$x[,1]),PC2=as.numeric(pca$x[,2]),tissue=color.var[,"TISSUE"],strain=color.var[,"STRAIN"],age=color.var[,"AGE"],gender=color.var[,"GENDER"])
#p=vector("list",4)
for(i in 3:ncol(dat))
  {

p <- ggplot(dat, aes(PC1,PC2))+ geom_point(aes(color=dat[,i]))+labs(x=xl,y=yl,color=colnames(dat)[i])
print(p)

  }
#multiplot(p[[1]],p[[2]],p[[3]],p[[4]],cols=2) ## this doesn't work
}

regression <- function(dat)
  {
    percent.vs.age=matrix(NA,nr=ncol(dat)-5,nc=4)
    for(i in 6:ncol(dat))
      {
        temp=summary(lm(dat[,i]~dat$age))$coef
        if(nrow(temp)==2) percent.vs.age[i-5,]=temp[,c(1,4)]
      }

    percent.vs.age= percent.vs.age[,-ncol(percent.vs.age)+1]

    rownames(percent.vs.age)=colnames(dat)[6:ncol(dat)]
    colnames(percent.vs.age)=c("intercept","aging_coef","aging_p")

    return(percent.vs.age)
  }

change.column.name <- function(x)
  {
    ## y=x
    ## for(i in 1:length(x))
    ##   {
    ##     temp= strsplit(x[i],"[.]")[[1]]
    ##     temp[temp==
    ##     y[i]=paste(temp[-c(1:6,length(temp)-0:3)],collapse=" ")
    ##   }
    return(x)
  }
pvalue.convert <- function(x)
  {
    for(i in 1:ncol(x))
      {
        x[,i]=p.adjust(x[,i],"fdr")
      }
    y=x
    y[]=0
    y[x<0.05]=1
    return(y)
    
  }

    
####### CD4T,8T subpopulation ggplot ##############



plotting <- function(dat0,datF,datM,datB6,datNZO)
  {
    for(k in 6:ncol(dat0))
      { 
        rg=range(dat0[,k],na.rm=T)
        rgx=range(dat0$age,na.rm=T)

        pF <- ggplot(datF, aes(age, datF[,k] ))+geom_point(aes(color=type)) +ggtitle(paste(colnames(dat)[k],"Female"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        pM <- ggplot(datM, aes(age, datM[,k] ))+geom_point(aes(color=type)) +ggtitle(paste(colnames(dat)[k],"Male"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

        pB6 <- ggplot(datB6, aes(age, datB6[,k] ))+geom_point(aes(color=Sex)) +ggtitle(paste(colnames(dat)[k],"B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        pNZO <- ggplot(datNZO, aes(age, datNZO[,k] ))+geom_point(aes(color=Sex)) +ggtitle(paste(colnames(dat)[k],"NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

        multiplot(pF,pM,pB6,pNZO,cols=2)
      }
  }

plotting2 <- function(dat0,datMB6,datMNZO,datFB6,datFNZO)
{
  for(k in 6:ncol(dat0))
    {
        rg=range(dat0[,k],na.rm=T)
        rgx=range(dat0$age,na.rm=T)

        pFB6 <- ggplot(datFB6, aes(age, datFB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"F B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        pMNZO <- ggplot(datMNZO, aes(age, datMNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"M NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

        pMB6 <- ggplot(datMB6, aes(age, datMB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"M B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        pFNZO <- ggplot(datFNZO, aes(age, datFNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"F NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)

        
        multiplot(pFB6,pFNZO,pMB6,pMNZO,cols=2)
      }
}

plotting3 <- function(dat0,datB6,datNZO)
{
  for(k in 6:ncol(dat0))
    {
        rg=range(dat0[,k],na.rm=T)
        rgx=range(dat0$age,na.rm=T)

        pB6 <- ggplot(datB6, aes(age, datB6[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"B6"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        pNZO <- ggplot(datNZO, aes(age, datNZO[,k] ))+geom_point() +ggtitle(paste(colnames(dat)[k],"NZO"))+ylab("percentage")+xlim(rgx[1],rgx[2])+ylim(rg[1],rg[2])+geom_smooth(method='lm',se=F)
        
        multiplot(pB6,NA,pNZO,NA,cols=2)
      }
}

subplot <- function(ggdf)
  {
    ggplot(ggdf, aes(x = new, y = value, fill = variable)) +
      geom_bar(stat = "identity") +  facet_grid(tissue~ type, scales = "free_x") + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Age - sample ID")+ylab("percentage") #   +  scale_fill_manual(values = c("red","blue","green"))
  }


subplot <- function(ggdf,title)
  {

    ggplot(ggdf, aes(x = new, y = value, colour = type)) +
      geom_jitter(width = 0.1, height = 0.1) +  facet_grid(tissue~ variable, scales = "free_x") + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("age")+ylab(paste("percentage"))+ggtitle(title)+geom_smooth(method = "lm", se = FALSE,size=0.3)+ scale_x_continuous(breaks = c(3,12,18))# stat_smooth(method="lm",col="grey", se=FALSE)#   +  scale_fill_manual(values = c("red","blue","green"))

  }


subplot2 <- function(ggdf,marker)
  {

    temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),marker)))
    temp=temp[nrow(temp),]
    ggdf$variable=temp
    
    ggplot(ggdf, aes(x = new, y = value, fill = variable)) +
      geom_point(aes(colour=variable)) +  facet_grid(tissue~ type, scales = "free_x") + theme_bw() 
      + scale_x_continuous(breaks = c(3,12,18))+  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("Age - sample ID")+ylab(paste("percentage of",marker,"+")) #   +  scale_fill_manual(values = c("red","blue","green"))
  }
  

subplot3 <- function(ggdf,marker)
  {
    temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),marker)))
    temp=temp[nrow(temp),]
    ggdf$variable=temp

    ggplot(ggdf, aes(x = new, y = value,colour=variable)) +
          geom_jitter(width = 0.1, height = 0.1) +  facet_grid(tissue~ type, scales = "free_x")+xlab("Age")+ylab("") + theme_bw() +        theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_color_manual(values = c("red","blue","green","yellow"),   breaks=c("X..Viable","Granulocytes.Fre.Viable.Cells", "No.Granulocytes.Monocytes.Fre.Parent", "No.Granulocytes.Monocytes.Fre.Viable.Cells"),   labels=c("Viable.Cells.percentage","Granulocytes.\nFreq.of.Viable.Cells", "No.Granulocytes.Monocytes.\nFreq.of.Parent", "No.Granulocytes.Monocytes.\nFreq.of.Viable.Cells")) #   +  scale_fill_manual(values = c("red","blue","green"))
  }
  
subplot4 <- function(ggdf,marker,title)
  {

    temp=as.matrix(data.frame(strsplit(as.character(ggdf$variable),paste(marker,".",sep=""))))
    temp=temp[nrow(temp),]
    ggdf$variable=temp
    
    ggplot(ggdf, aes(x = new, y = value, colour = type)) +
      geom_jitter(width = 0.1, height = 0.1) +  facet_grid(tissue~ variable, scales = "free_x") + theme_bw() +ggtitle(title)+geom_smooth(method = "lm", se = FALSE,size=0.3)+ scale_x_continuous(breaks = c(3,12,18))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("age")+ylab(paste("percentage of",marker,"+"))#   +  scale_fill_manual(values = c("red","blue","green"))
  }
