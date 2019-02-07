library(corrplot)

a=read.csv("../results/summaryATACResults.csv")
summary(a[,-(1:2)])
x=cor(a[,-(1:2)])

corrplot(x,type="upper",diag=FALSE,tl.srt=45)
#mtext("correlation", at=3, line=-8, cex=1)
