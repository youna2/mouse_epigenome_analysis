
### remove problematic sample ###

#load("../../flow/src/newproblemsampleID.Rdata")
selectsample= !(as.numeric(SAMPLEMOUSEID) %in% c(26,27,250,252))#,problemsample))
Y=Y[,selectsample]
bed=(bed[,selectsample])

### keep genes/peaks with cpm>1 for at least 2 samples
selectgene= rowSums(bed>1)>=2
bed=bed[selectgene,]
Y=Y[selectgene,]
annotation=cbind(annotation[selectgene,])

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
    selectgene= annotation[,"Nearest.Ensembl"] %in% exp.gene & abs(as.numeric(annotation[,"Distance.to.TSS"]))<1000

    
    filename="../results/ATACseqData2.Rdata"
  }


bed=log(bed+min(bed[bed>0]))
AGE=AGE[selectsample]
GENDER=GENDER[selectsample]
TISSUE=TISSUE[selectsample]
STRAIN=STRAIN[selectsample]
TYPE=TYPE[selectsample]
SAMPLEMOUSEID=SAMPLEMOUSEID[selectsample]

if(sum(strsplit(getwd(),"/")[[1]]=="RNAseq")==1 ) save(bed,Y,annotation,AGE,GENDER,TISSUE,STRAIN,TYPE,SAMPLEMOUSEID,file=filename)
