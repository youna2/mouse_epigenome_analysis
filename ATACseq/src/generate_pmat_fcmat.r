
p.cutoff=0.05

p.mat.strain=NULL ##bed[,"p.value"]
fc.mat.strain= NULL  ##-bed[,"Fold"]
if(tid>=5)
  {
    commonpattern=TRUE
    if(tid==5) int.or.slope="B6 vs NZO common"
    if(tid==6) int.or.slope="B6 vs NZO opposite"
    if(tid==7) int.or.slope="only NZO significant"
    if(tid==8) int.or.slope="only B6 significant"
    if(tid==9) int.or.slope="within B6"
    if(tid==10) int.or.slope="within NZO"
    if(tid==11) int.or.slope="all"
    
  }else{
    commonpattern=FALSE
    if(tid==3)
      {
        typeinteraction=TRUE
        int.or.slope="slope"
      }else{
        typeinteraction=FALSE
        int.or.slope="intercept"
      }
  }
#### remove samples of age 26.75 months ######

selectsample= (AGE != 26.75)
y0=Y[,selectsample]
y0forheatmap=bed[,selectsample]
age0=AGE[selectsample]
gender0=GENDER[selectsample]
tissue0=TISSUE[selectsample]
type0=TYPE[selectsample]
samplemouseid0=SAMPLEMOUSEID[selectsample]

## Do differential analysis between NZO and B6 in each tissue type ######

setwd("../results")
previous.dir=getwd()
utissue=c( "spleen", "BM", "memory","naive", "PBL")         

if(commonpattern)
  {
    pdf(paste(int.or.slope,"pattern.pdf",sep="_"))
    topgene=int.or.slope
    if(tid==5 | tid>=9)
      {
        positivepeak="opening"
        negativepeak="closing"
      }else{
        if(tid==8)
          {
            positivepeak="higher in B6"
            negativepeak="lower in B6"
          }else{
            positivepeak="higher in NZO"
            negativepeak="lower in NZO"
          }
      }
  }else{

    positivepeak="higher in NZO"
    negativepeak="lower in NZO"

    if(typeinteraction)
      {
        pdf("diffpeakstrain_B6_vs_NZO_slope.pdf")
        topgene="B6_vs_NZO_slope"
      }else{
        pdf("diffpeakstrain_B6_vs_NZO_intercept.pdf")
        topgene="B6_vs_NZO_intercept"
      }
  }





for(i in 1:length(utissue))
  {
    y=y0[,tissue0==utissue[i]]
    age=age0[tissue0==utissue[i]]
    gender=gender0[tissue0==utissue[i]]
    type=type0[tissue0==utissue[i]]
    if(length(unique(age[type=="B6"]))>1 & length(unique(age[type=="NZO"]))>1)
      {
### Do differential or common analysis between NZO and B6

        
        if(commonpattern)  atac.glmtop.strain=edgeRfitcommon(y,age,gender,type) else atac.glmtop.strain=edgeRfitstrain(y,age,gender,type)
        

        sel=which(atac.glmtop.strain[,"FDR"]<p.cutoff)
#        write.table(cbind(annotation[ sel,"Entrez.ID"],NA,atac.glmtop.strain[sel,"logFC"]),file=paste(topgene,paste(utissue[i],int.or.slope),".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

        
        
#### draw heatmap of age-increasing/decreasing peaks or genes  #################
        if(length(sel)>5)
          {
            atac.glmtop=atac.glmtop.strain
            p.mat.strain=cbind(p.mat.strain,atac.glmtop[,"FDR"])
            fc.mat.strain=cbind(fc.mat.strain,atac.glmtop[,"logFC"])
            heatmapmat=y0forheatmap[atac.glmtop[,"FDR"]<p.cutoff,tissue0==utissue[i]]

            
            annot.row=annot=atac.glmtop[atac.glmtop[,"FDR"]<p.cutoff,"logFC"]

            annot[annot.row>0]=positivepeak
            annot[annot.row<0]=negativepeak
            
            
            
            colnames(heatmapmat)=paste( type,"_", samplemouseid0[tissue0==utissue[i]] ,"_",age,sep="")
            

            heatmapmat=heatmapmat[,order(age)]
            type=type[order(age)]
            heatmapmat=heatmapmat[,order(type)]


            if(nrow(heatmapmat)>40000)
              {
                tmpsample=sample(1:nrow(heatmapmat),40000)
                heatmapmat=heatmapmat[tmpsample,]
                annot=annot[tmpsample]
              }

            annot=as.data.frame(annot)
            rownames(annot)=rownames(heatmapmat)

            if(commonpattern)
              {
                colnames(annot)="common"
              }else{        
                if(typeinteraction)
                  {
                    colnames(annot)="slope"
                  }else{
                    colnames(annot)="intercept"
                  }
              }
            if(min(heatmapmat)==0) heatmapmat=heatmapmat+1
            
           # pheatmap(log(heatmapmat),scale="row",cluster_cols = FALSE,main=paste(utissue[i], sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]>0),positivepeak, sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]<0),negativepeak),annotation_row=annot,show_rownames=F,color=colorRampPalette(c("blue","white","red"))(100))  
                                        #   
            colnames(fc.mat.strain)[ncol(fc.mat.strain)]=colnames(p.mat.strain)[ncol(p.mat.strain)]=paste(utissue[i],int.or.slope)
          }
      }
  }
save(p.mat.strain,fc.mat.strain,file=paste("pmat_fcmat_",int.or.slope,".Rdata",sep=""))

dev.off()
