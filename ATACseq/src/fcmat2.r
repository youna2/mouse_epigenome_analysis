source("../../ATACseq/src/library_function.r")


load("../../ATACseq/results3/ATACseqData2.Rdata")
load("../../ATACseq/results3/pmat_fcmat_within B6.Rdata")
p.mat.B6=p.mat.strain
fc.mat.B6=fc.mat.strain
meanexp.B6=meanexp

load("../../ATACseq/results3/pmat_fcmat_within NZO.Rdata")
p.mat.NZO=p.mat.strain
fc.mat.NZO=fc.mat.strain
meanexp.NZO=meanexp

humangenename=rownames(annotation)

save(meanexp.B6,fc.mat.B6,p.mat.B6,meanexp.NZO,fc.mat.NZO,p.mat.NZO,humangenename,file="../../ATACseq/results/mouseATACFC.Rdata")   
