##This is to compare FC as well as mean expression of inflammation genes between B6 and NZO##
pdf(file="temp_inflammation_mean.pdf")
load("/Users/youna/mouse_aging_data_analysis_code/ATACseq/results/inflammation1_gene_list.Rdata")
inflamgene=toupper(convert.mouse.entrez.to.human.symbol(rownames(bed))) %in% toupper(inflammation1)
bed.inflammation=bed[inflamgene,]
utissue=unique(TISSUE)
uage=c(3,18)
par(mfrow=c(2,2))
for(i in 1:length(utissue))
  {
    for(k in 1:length(uage))
      {### compare mean B6 vs mean NZO 
        meanB6=rowMeans(bed.inflammation[,AGE==uage[k]&TISSUE==utissue[i]&STRAIN=="B6"])
        meanNZO=rowMeans(bed.inflammation[,AGE==uage[k]&TISSUE==utissue[i]&STRAIN=="NZO"])
        plot(meanB6,meanNZO,main=paste(utissue[i],uage[k],"months"));abline(0,1)
      }
  }
load("../results/positive_inflammation_gene.Rdata")
bed.inflammation=bed.inflammation[rownames(bed.inflammation)%in%positive.inflammation,]
rownames(bed.inflammation)=convert.mouse.entrez.to.human.symbol(rownames(bed.inflammation))
##


library(ComplexHeatmap)                                

### This is to compare FC of B6 vs NZO for RNA and ATAC
par(mfrow=c(2,2))
if(sum(strsplit(getwd(),"/")[[1]]=="RNAseq")==1)
  {
    load("/Users/youna/mouse_aging_data_analysis_code/RNAseq/results/pmat_fcmat_within NZO_before_filtering.Rdata")
    rownames(fc.mat.strain)=annotation.orig[,"Entrez.ID"]
    inflamgene=toupper(convert.mouse.entrez.to.human.symbol(rownames(fc.mat.strain))) %in% toupper(inflammation1)
    fc.NZO=fc.mat.strain[inflamgene,]
    p.NZO=p.mat.strain[inflamgene,]
    load("/Users/youna/mouse_aging_data_analysis_code/RNAseq/results/pmat_fcmat_within B6_before_filtering.Rdata")
    fc.B6=fc.mat.strain[inflamgene,]
    p.B6=p.mat.strain[inflamgene,]

    for(i in 1:4)
      {plot(fc.NZO[,i],fc.B6[,i],main=strsplit(colnames(fc.mat.strain)[i]," ")[[1]][1],xlab="FC NZO",ylab="FC B6");abline(0,1)}

    positive.inflammation=rownames(fc.NZO)[rowSums(fc.NZO[,1:2]>0.01)>=2 & rowSums(fc.B6[,1:2]>0.01)>=2]
    save(positive.inflammation,file="../results/positive_inflammation_gene.Rdata")
  }else{

    positive.inflammation=vector("list",5)
    
    for(BTID in 2:5)
      {
        load(paste("/Users/youna/mouse_aging_data_analysis_code/ATACseq/results",BTID,"/pmat_fcmat_within NZO_before_filtering.Rdata",sep=""))
        rownames(fc.mat.strain)=annotation.orig[,"Entrez.ID"]
        inflamgene=toupper(convert.mouse.entrez.to.human.symbol(rownames(fc.mat.strain))) %in% toupper(inflammation1)
        fc.NZO=cbind(fc.mat.strain[inflamgene,])
        p.NZO=cbind(p.mat.strain[inflamgene,])

        load(paste("/Users/youna/mouse_aging_data_analysis_code/ATACseq/results",BTID,"/pmat_fcmat_within B6_before_filtering.Rdata",sep=""))
        fc.B6=cbind(fc.mat.strain[inflamgene,])
        p.B6=cbind(p.mat.strain[inflamgene,])

        plot(fc.NZO,fc.B6,main=strsplit(colnames(fc.mat.strain)[1]," ")[[1]][1],xlab="FC NZO",ylab="FC B6");abline(0,1)

        positive.inflammation[[BTID]]=rownames(fc.NZO)[fc.NZO>0.02 & fc.B6>0.02]
      }

    positive.inflammation=intersect(positive.inflammation[[2]],positive.inflammation[[3]])

    save(positive.inflammation,file="../results/positive_inflammation_gene.Rdata")

  }
take.mean <- function(x)
  {
    gene=unique(rownames(x))
    y=matrix(NA,nr=length(gene),nc=ncol(x))
    rownames(y)=gene
    for(i in 1:length(gene))
      {
        y[i,]=colMeans(rbind(x[rownames(x)==gene[i],]))
      }
    return(y)
  }
### This is to draw heatmap of gene or peak level of genes in inflammation module           
bed.inflammation=take.mean(bed.inflammation)
bed.inflammation=bed.inflammation[,order(STRAIN,AGE)]
AGEi=AGE[order(STRAIN,AGE)]
STRAINi=STRAIN[order(STRAIN,AGE)]
TISSUEi=TISSUE[order(STRAIN,AGE)]

for(i in 1:length(utissue))
draw(Heatmap(bed.inflammation[,TISSUEi==utissue[i]],cluster_columns=F,heatmap_legend_param = list(title = "log"),column_title=utissue[i], row_names_gp = gpar(fontsize = 9),row_title="",show_column_names=F,top_annotation = HeatmapAnnotation(age = AGEi[TISSUEi==utissue[i]],strain=STRAINi[TISSUEi==utissue[i]],   col = list(age=c("3"="blue","12"="white","18"="red"),strain = c("B6" =  "green", "NZO" = "yellow"))   )))

dev.off()
####
