
p.cutoff=0.1

meanexp=p.mat.strain=NULL ##bed[,"p.value"]
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

y0=Y

age0=AGE
gender0=GENDER
tissue0=TISSUE
type0=TYPE
samplemouseid0=SAMPLEMOUSEID
librarysize0=LIBRARYSIZE
## Do differential analysis between NZO and B6 in each tissue type ######

setwd(paste("../results",BTID,sep=""))
previous.dir=getwd()
utissue=unique(tissue0)         

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
#    if(tid==9) meanexp=cbind(meanexp,rowMeans(bed[,tissue0==utissue[i] & type0=="B6"]))
    
#    if(tid==10) meanexp=cbind(meanexp,rowMeans(bed[,tissue0==utissue[i] & type0=="NZO"]))
    
    y=y0[,tissue0==utissue[i]]
    age=age0[tissue0==utissue[i]]
    gender=gender0[tissue0==utissue[i]]
    type=type0[tissue0==utissue[i]]
    librarysize=librarysize0[tissue0==utissue[i]]

    
#    if(length(unique(age[type=="B6"]))>1 & length(unique(age[type=="NZO"]))>1)
    if(TRUE)
      {
### Do differential or common analysis between NZO and B6

        
        if(commonpattern)  atac.glmtop.strain=edgeRfitcommon(y,age,gender,type,librarysize) else atac.glmtop.strain=edgeRfitstrain(y,age,gender,type)

        if(tid==9 | tid==10) MAplot(atac.glmtop.strain,utissue[i])

        sel=which(atac.glmtop.strain[,"FDR"]<p.cutoff)
#        write.table(cbind(annotation[ sel,"Entrez.ID"],NA,atac.glmtop.strain[sel,"logFC"]),file=paste(topgene,paste(utissue[i],int.or.slope),".txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")



        meanexp=cbind(meanexp,atac.glmtop.strain[,"logCPM"])
        


        
#### draw heatmap of age-increasing/decreasing peaks or genes  #################
            atac.glmtop=atac.glmtop.strain

        if(sum(strsplit(getwd(),"/")[[1]]=="ATACseq")==1 )
          {

        if(tid==9)
          {
            write.table(annotation.orig[atac.glmtop[,"FDR"]<0.05 & atac.glmtop[,"logFC"]>0  ,1:3],file=paste("diffpeak_pos_B6_",utissue[i],".bed",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
            write.table(annotation.orig[atac.glmtop[,"FDR"]<0.05 & atac.glmtop[,"logFC"]<0  ,1:3],file=paste("diffpeak_neg_B6_",utissue[i],".bed",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
          }  


        if(tid==10)
          {
            write.table(annotation.orig[atac.glmtop[,"FDR"]<0.05 & atac.glmtop[,"logFC"]>0  ,1:3],file=paste("diffpeak_pos_NZO_",utissue[i],".bed",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
            write.table(annotation.orig[atac.glmtop[,"FDR"]<0.05 & atac.glmtop[,"logFC"]<0  ,1:3],file=paste("diffpeak_neg_NZO_",utissue[i],".bed",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
          }  
      }

        if(tid==9) write.table(cbind(annotation.orig,atac.glmtop),file=paste("supplementary_B6_FDR_FC_",utissue[i],".txt",sep=""),quote=F,row.names=F,sep="\t")
        if(tid==10) write.table(cbind(annotation.orig,atac.glmtop),file=paste("supplementary_NZO_FDR_FC_",utissue[i],".txt",sep=""),quote=F,row.names=F,sep="\t")
        p.mat.strain=cbind(p.mat.strain,atac.glmtop[,"FDR"])
        fc.mat.strain=cbind(fc.mat.strain,atac.glmtop[,"logFC"])

        if(length(sel)>5)
          {

            heatmapmat=y0forheatmap[atac.glmtop[,"FDR"]<p.cutoff,tissue0==utissue[i]]

            
            annot.row=annot=atac.glmtop[atac.glmtop[,"FDR"]<p.cutoff,"logFC"]

            annot[annot.row>0]=positivepeak
            annot[annot.row<0]=negativepeak
            
            
            
            colnames(heatmapmat)=paste( type,"_", samplemouseid0[tissue0==utissue[i]] ,"_",age,sep="")
            

            heatmapmat=heatmapmat[,order(age)]
            type=type[order(age)]
            heatmapmat=heatmapmat[,order(type)]


            if(nrow(heatmapmat)>20000)
              {
                tmpsample=sample(1:nrow(heatmapmat),20000)
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
            if(min(heatmapmat)==0) heatmapmat=heatmapmat+min(heatmapmat[heatmapmat>0])
            
            pheatmap(log(heatmapmat),scale="none",cluster_cols = FALSE,main=paste(utissue[i], sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]>0),positivepeak, sum(atac.glmtop[,"FDR"]<p.cutoff & atac.glmtop[,"logFC"]<0),negativepeak),annotation_row=annot,show_rownames=F,color=colorRampPalette(c("blue","white","red"))(100))  
                                        #   
      }
        colnames(fc.mat.strain)[ncol(fc.mat.strain)]=colnames(p.mat.strain)[ncol(p.mat.strain)]=paste(utissue[i],int.or.slope)
          
      }
  }
save(atac.glmtop,annotation.orig,meanexp,p.mat.strain,fc.mat.strain,file=paste("pmat_fcmat_",int.or.slope,"_before_filtering.Rdata",sep=""))

print("dim bed, Y, annotation, pmat before")
print(dim(bed))
print(dim(Y))
print(dim(annotation))
print(dim(p.mat.strain))
if(sum(strsplit(getwd(),"/")[[1]]=="ATACseq")==1 )
  {
    tmp=colnames(p.mat.strain)
    p.mat.strain=cbind(p.mat.strain[selectgene,])
    colnames(p.mat.strain)=tmp
    tmp=colnames(fc.mat.strain)
    fc.mat.strain=cbind(fc.mat.strain[selectgene,])
    colnames(fc.mat.strain)=tmp

    if(tid==9 | tid==10) meanexp=cbind(meanexp[selectgene,])
  }

save(annotation,meanexp,p.mat.strain,fc.mat.strain,file=paste("pmat_fcmat_",int.or.slope,".Rdata",sep=""))

dev.off()
print("dim bed, Y, annotation, pmat after (only ATAC change)")
print(dim(bed))
print(dim(Y))
print(dim(annotation))
print(dim(p.mat.strain))
