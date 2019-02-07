med=quantile(as.vector(bed),0.95)
quantile(rowMeans(bed>med),0:10*0.1)
sum(rowMeans(bed>med)>0.95)
sel= rowMeans(bed>med)>0.95
sort(colMeans(bed[sel,]>med))
