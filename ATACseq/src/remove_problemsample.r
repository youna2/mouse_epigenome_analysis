
### remove problematic sample ###

load("../../flow/src/newproblemsampleID.Rdata")

table(SAMPLEMOUSEID[SAMPLEMOUSEID%in%problemsample])
selectsample= !(  (as.numeric(SAMPLEMOUSEID) %in% c(26,27,250,252,problemsample)) | (SAMPLEMOUSEID==36&TISSUE=="CD8.memory")  | (SAMPLEMOUSEID==36&TISSUE=="CD8.naive") | (SAMPLEMOUSEID==39&TISSUE=="CD8.naive") | (SAMPLEMOUSEID==42&TISSUE=="CD8.memory"&STRAIN=="NZO") )

#selectsample= !(as.numeric(SAMPLEMOUSEID) %in% c(26,27,250,252))#,problemsample))
Y=Y[,selectsample]
bed=(bed[,selectsample])

### keep genes/peaks with cpm>1 for at least 2 samples and for ATACseq data peaks whose peak width is larger than the 20 percentile
if(sum(strsplit(getwd(),"/")[[1]]=="ATACseq")==1 )
{
peakwidth= as.numeric(annotation[,"End"])-as.numeric(annotation[,"Start"])
selectgene= rowSums(bed>1)>=2 & peakwidth > quantile(peakwidth,0.2) 
}else selectgene= rowSums(bed>1)>=2

bed=bed[selectgene,]
Y=Y[selectgene,]

annotation=cbind(annotation[selectgene,])
annotation.orig=annotation
### save these gene list

if(sum(strsplit(getwd(),"/")[[1]]=="RNAseq")==1 )
  {
    colnames(annotation)="Entrez.ID"
    exp.gene=rownames(bed)
    save(exp.gene,file="../results/expressed_genelist.Rdata")
    filename="../results/RNAseqData2.Rdata"
  }

### Only keep peaks that are annotated to the expressed genes

if(sum(strsplit(getwd(),"/")[[1]]=="ATACseq")==1 )
  {    
    load("../../RNAseq/results/expressed_genelist.Rdata")
    selectgene= annotation[,"Nearest.Ensembl"] %in% exp.gene & abs(as.numeric(annotation[,"Distance.to.TSS"]))< DistanceToTSS

    annotation=annotation[selectgene,]
       
    filename=paste("../results",BTID,"/ATACseqData2.Rdata",sep="")
  }

y0forheatmap=bed
bed=log(bed+min(bed[bed>0]))
AGE=AGE[selectsample]
GENDER=GENDER[selectsample]
TISSUE=TISSUE[selectsample]
STRAIN=STRAIN[selectsample]
TYPE=TYPE[selectsample]
SAMPLEMOUSEID=SAMPLEMOUSEID[selectsample]
LIBRARYSIZE=LIBRARYSIZE[selectsample]
save(LIBRARYSIZE,annotation.orig,annotation,AGE,GENDER,TISSUE,STRAIN,TYPE,SAMPLEMOUSEID,file=filename)
