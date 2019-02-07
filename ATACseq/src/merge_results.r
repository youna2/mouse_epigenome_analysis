number=vector("list",10)
barplot.all=balloonplot.all=vector("list",4)
tissue=c("BM","spleen","PBL","memory","naive")

for(BTID in 2:5)
  {
    print(tissue[BTID])
    load(paste("../../ATACseq/data/ATACseqData",BTID,".Rdata",sep=""))
    print(mean(as.numeric(bed[,3])+1==as.numeric(annotation[,2])))
    print(mean(as.numeric(bed[,4])==as.numeric(annotation[,3])))
    print(mean(bed[,2]==annotation[,1]))
    print(dim(bed))

    
    load(paste("../../ATACseq/results",BTID,"/balloonplot",BTID,".Rdata",sep=""))
    for(i in 5:10)
      {
        if(!is.null(total.number.plot[[i]][[1]]))
          number[[i]]=rbind(number[[i]],total.number.plot[[i]][[1]])
      }
    
    for(i in 1:4)
      {
        if(!is.null(balloonplot[[i]]))
          balloonplot.all[[i]]=rbind(balloonplot.all[[i]],balloonplot[[i]])
        if(!is.null(barplot[[i]]))
          barplot.all[[i]]=rbind(barplot.all[[i]],barplot[[i]][barplot[[i]][,5]!="CD8.memory" & barplot[[i]][,5]!="CD8.naive",])

      }
  }

balloonplot=balloonplot.all
barplot=barplot.all
pdf(file="../../ATACseq/results/summary_pathwayplot.pdf",height=17,width=10)
plot.balloonplot.mouse(balloonplot)    
draw.barplot(barplot)
diff.peaks()
dev.off()

save(balloonplot,barplot,file="../../ATACseq/results/balloonplotATAC.Rdata")
ylimmax=15000
YLIM=c(-ylimmax,ylimmax)
int.or.slope=c("","","slope","intercept","B6 vs NZO common","B6 vs NZO opposite","only NZO significant","only B6 significant","within B6","within NZO")


  
pdf(file="../../ATACseq/results/summaryfigure.pdf")
temp=vector("list",length(number))
for(i in 5:10)
  temp[[i]]=twoway.barplot(number[[i]][,1],number[[i]][,2],number[[i]][,3],(-10):10*1000,(-10):10*1000,"Tissue","no. differential peaks/genes",int.or.slope[i],YLIM)[[2]]


multiplot(temp[[9]],temp[[5]],temp[[8]],temp[[4]],temp[[10]],temp[[6]],temp[[7]],temp[[3]],cols=2)
dev.off()
