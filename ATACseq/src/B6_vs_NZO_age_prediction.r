    title="NZO"
        sel=which(STRAIN==title)
        x_training=t(bed)[-sel,]
        y_training=AGE[-sel]

        x_testing=rbind(t(bed)[sel,])
        y_testing=AGE[sel]


        fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
        fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
        fit_ridge <- cv.glmnet(x=x_training, y=y_training,  alpha=0)


        pred_lasso <- predict(fit_lasso, s="lambda.min", newx=x_testing)
        pred_elastic <- predict(fit_elastic, s="lambda.min", newx=x_testing)
        pred_ridge <- predict(fit_ridge, s="lambda.min", newx=x_testing)
pdf(file="B6_vs_NZO_age_prediction.pdf")
par(mfrow=c(3,3))
plot(pred_lasso,y_testing,main=paste("testing =",title))
abline(a=0,b=1)
plot(pred_elastic,y_testing,main=paste("testing =",title))
abline(a=0,b=1)
plot(pred_ridge,y_testing,main=paste("testing =",title))
abline(a=0,b=1)

title="B6"

        sel=which(STRAIN==title)
        x_training=t(bed)[-sel,]
        y_training=AGE[-sel]

        x_testing=rbind(t(bed)[sel,])
        y_testing=AGE[sel]


        fit_lasso <- cv.glmnet(x=x_training, y=y_training,  alpha=1)
        fit_elastic <- cv.glmnet(x=x_training, y=y_training, alpha=0.5)
        fit_ridge <- cv.glmnet(x=x_training, y=y_training,  alpha=0)


        pred_lasso <- predict(fit_lasso, s="lambda.min", newx=x_testing)
        pred_elastic <- predict(fit_elastic, s="lambda.min", newx=x_testing)
        pred_ridge <- predict(fit_ridge, s="lambda.min", newx=x_testing)


plot(pred_lasso,y_testing,main=paste("testing =",title))
abline(a=0,b=1)
plot(pred_elastic,y_testing,main=paste("testing =",title))
abline(a=0,b=1)
plot(pred_ridge,y_testing,main=paste("testing =",title))
abline(a=0,b=1)


dev.off()
