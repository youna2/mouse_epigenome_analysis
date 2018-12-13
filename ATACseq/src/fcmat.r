source("../../ATACseq/src/library_function.r")


load("../../ATACseq/results/ATACseqData2.Rdata")
load("../../ATACseq/results/pmat_fcmat_within B6.Rdata")
p.mat.B6=p.mat.strain
fc.mat.B6=fc.mat.strain

load("../../ATACseq/results/pmat_fcmat_within NZO.Rdata")
p.mat.NZO=p.mat.strain
fc.mat.NZO=fc.mat.strain


onlypromoter=T
if(onlypromoter)
  {
#    sel=(removing.paren(annotation[,"Annotation"])=="promoter-TSS ")
    sel=(abs(as.numeric(annotation[,"Distance.to.TSS"]))<1000)
    annotation=annotation[sel,]
    fc.mat.B6=fc.mat.B6[sel,]
    p.mat.B6=p.mat.B6[sel,]
    fc.mat.NZO=fc.mat.NZO[sel,]
    p.mat.NZO=p.mat.NZO[sel,]
  }


###

        EntrezATAC=annotation[,"Entrez.ID"]

        common=unique(EntrezATAC)
        length(common)
        matching.index=matrix(NA,nr=length(common),nc=1)
        argmin <- function(x) return(which(x==min(x))[1])

        for(j in 1:length(common))
          {
            temp=which(EntrezATAC %in% common[j])
            x=argmin(abs(as.numeric(annotation[temp,"Distance.to.TSS"])))
            matching.index[j,1]=temp[x]
          }
humangenename=convert.mouse.entrez.to.human.symbol(annotation[matching.index[,1],"Entrez.ID"])

fc.mat.B6=fc.mat.B6[matching.index[,1],]
p.mat.B6=p.mat.B6[matching.index[,1],]
fc.mat.NZO=fc.mat.NZO[matching.index[,1],]
p.mat.NZO=p.mat.NZO[matching.index[,1],]
###



save(fc.mat.B6,p.mat.B6,fc.mat.NZO,p.mat.NZO,humangenename,file="../../ATACseq/results/mouseATACFCpromoter.Rdata")   

###########################


load("../../RNAseq/results/RNAseqData2.Rdata")
load("../../RNAseq/results/pmat_fcmat_within B6.Rdata")
fc.mat.B6=fc.mat.strain
p.mat.B6=p.mat.strain

load("../../RNAseq/results/pmat_fcmat_within NZO.Rdata")
fc.mat.NZO=fc.mat.strain
p.mat.NZO=p.mat.strain




humangenename=convert.mouse.entrez.to.human.symbol(annotation[,1])


save(fc.mat.B6,p.mat.B6,fc.mat.NZO,p.mat.NZO,humangenename,file="../../RNAseq/results/mouseRNAFC.Rdata")   
     
## comparing predictors with the differential genes/peaks ###
x_training=t(bed)
rownames(x_training)=1:nrow(x_training)
y_training=AGE
fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
coef.est=coef(fit_elastic)[-1,]
length(coef.est)
dim(bed)
sum(coef.est!=0)

load("../../ATACseq/results/pmat_fcmat_within B6.Rdata")
p.mat.B6=p.mat.strain
fc.mat.B6=fc.mat.strain
 
load("../../ATACseq/results/pmat_fcmat_within NZO.Rdata")
p.mat.NZO=p.mat.strain
fc.mat.NZO=fc.mat.strain
table(rowSums(fc.mat.B6[coef.est>0,]>0)+rowSums(fc.mat.NZO[coef.est>0,]>0))
table(rowSums(fc.mat.B6[coef.est<0,]<0)+rowSums(fc.mat.NZO[coef.est<0,]<0))


#####
