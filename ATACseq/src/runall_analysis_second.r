
total.number.plot=bmpathimmune=spleenpathimmune=pblpathimmune=memorypathimmune= naivepathimmune=bmpathwiki=spleenpathwiki=pblpathwiki= memorypathwiki= naivepathwiki=bmpathkegg=spleenpathkegg=pblpathkegg= memorypathkegg= naivepathkegg=vector("list",10)

balloonplot=vector("list",4)
for(tid in 3:10)
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

pdf(file="summaryfigure.pdf")
temp=total.number.plot
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)

temp=bmpathimmune
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
temp=spleenpathimmune
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
temp=pblpathimmune
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
temp=memorypathimmune
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
temp= naivepathimmune
multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)


temp= bmpathwiki
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= spleenpathwiki
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= pblpathwiki
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= memorypathwiki
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= naivepathwiki
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])


temp= bmpathkegg
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= spleenpathkegg
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= pblpathkegg
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= memorypathkegg
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
temp= naivepathkegg
multiplot2(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]])
dev.off()
pdf(file="summary_pathwayplot.pdf",height=17,width=10)
for(i in 1:3)
  {
    multiballoonplot=vector("list",4)
    
    for(j in 1:4)
      {
        pathwaymat=balloonplot[[i]]
        if(i==3)
          {
            pathwaymat=pathwaymat[as.numeric(pathwaymat[,2])<10^(-5),]
            pathwaymat=pathwaymat[pathwaymat[,4]=="memory" | pathwaymat[,4]=="naive" ,]
          }
        
        if(i==1) pathwaymat=pathwaymat[pathwaymat[,4]!="memory" &pathwaymat[,4]!="naive" ,]
        
        if(j ==1) sel=pathwaymat[,5]=="within NZO" | pathwaymat[,5]=="within B6"
        if(i>2)
          {
            if(j ==2 | j==3) sel=pathwaymat[,5]=="B6 vs NZO common" | pathwaymat[,5]=="B6 vs NZO opposite"| pathwaymat[,5]=="only NZO significant" | pathwaymat[,5]=="only B6 significant"
          }else{
            if(j ==2) sel=pathwaymat[,5]=="B6 vs NZO common" | pathwaymat[,5]=="B6 vs NZO opposite"
            if(j==3) sel=pathwaymat[,5]=="only NZO significant" | pathwaymat[,5]=="only B6 significant"

          }
        if(j ==4) sel=pathwaymat[,5]=="intercept" | pathwaymat[,5]=="slope"        

        pathwaymat=pathwaymat[sel,]
    
        df=data.frame(pathid=pathwaymat[,1],logp=-log10(as.numeric(pathwaymat[,2])+10^(-30)),direction=(pathwaymat[,3]),tissue=pathwaymat[,4],class=pathwaymat[,5])

        multiballoonplot[[j]]=(ggballoonplot(df, x = "tissue", y = "pathid", size = "logp", fill = "direction", facet.by =c(".", "class") ,ggtheme = theme_bw()))+scale_fill_manual(values = c("down"="dodgerblue","up"="firebrick1"),guide=guide_legend(title = NULL)) #+
#  scale_fill_viridis_c(option = "C")
      }

    if(i<3)
      {
        multiplot(multiballoonplot[[1]],multiballoonplot[[2]],multiballoonplot[[3]],multiballoonplot[[4]],col=1)
      }else{
        multiplot(multiballoonplot[[1]],multiballoonplot[[2]],col=1);
        print(multiballoonplot[[4]]);

      }
  }
    
dev.off()
