
total.number.plot=bmpathimmune=spleenpathimmune=pblpathimmune=memorypathimmune= naivepathimmune=bmpathwiki=spleenpathwiki=pblpathwiki= memorypathwiki= naivepathwiki=bmpathkegg=spleenpathkegg=pblpathkegg= memorypathkegg= naivepathkegg=vector("list",10)

balloonplot=barplot=vector("list",4)
for(tid in 5:10)
  {

if(tid==1) selB6=TRUE
if(tid==2) selB6=FALSE


if(tid>2)
  {
    source("../../ATACseq/src/generate_pmat_fcmat.r")
    source("../../ATACseq/src/find_differential_peak_or_gene_strain.r")
 }else source("../../ATACseq/src/find_differential_peak_or_gene.r")
}


multiplot2 <- function(a,b,c,d,e,f,g,h)
  {
    multiplot(a,e,cols=1)
    multiplot(b,f,cols=1)
    multiplot(c,g,cols=1)
    multiplot(d,h,cols=1)    
  }

## pdf(file="summaryfigure.pdf")
## temp=total.number.plot
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)

## temp=bmpathimmune
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
## temp=spleenpathimmune
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
## temp=pblpathimmune
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
## temp=memorypathimmune
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
## temp= naivepathimmune
## multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)


## temp= bmpathwiki
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= spleenpathwiki
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= pblpathwiki
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= memorypathwiki
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= naivepathwiki
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])


## temp= bmpathkegg
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= spleenpathkegg
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= pblpathkegg
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= memorypathkegg
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## temp= naivepathkegg
## multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
## dev.off()
## pdf(file="summary_pathwayplot.pdf",height=17,width=10)
## plot.balloonplot.mouse(balloonplot)    
## dev.off()

### Here, only total.number.plot and balloonplot matter. Others are extra
save(total.number.plot,bmpathimmune,spleenpathimmune,pblpathimmune,memorypathimmune, naivepathimmune,bmpathwiki,spleenpathwiki,pblpathwiki, memorypathwiki, naivepathwiki,bmpathkegg,spleenpathkegg,pblpathkegg, memorypathkegg, naivepathkegg,balloonplot,barplot,file=paste("balloonplot",BTID,".Rdata",sep=""))
