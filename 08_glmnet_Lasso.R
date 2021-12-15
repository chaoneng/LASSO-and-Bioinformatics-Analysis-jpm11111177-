##' 2021-04-26
##' Training a logistic regression classification model
##' http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
##' https://rpubs.com/skydome20/R-Note18-Subsets_Shrinkage_Methods
##' https://cosx.org/2016/10/data-mining-1-lasso/

rm(list = ls())
Sys.setenv(LANGUAGE="en")
setwd("C:/Users/flora/Desktop/gynecological cancers")

## Loading required packages
library("glmnet")
library("dplyr")

#-----------------------------------------------------------------#
# 1. Preparing the data
#-----------------------------------------------------------------#

mydata <- read.csv("./02_table/exprSet_CC_EC_after_combat.txt", header = TRUE, sep = "\t", row.names = 1) #13228 35
mypdata <- read.csv("./02_table/Clinical_information_CC_EC.csv", header = TRUE, sep = ",", row.names = 1)

degs <- read.csv("./02_table/78_common_DE_genes.csv", row.names = 1)
degs <- rownames(degs) #78 differentially expressed genes


mydata <- t(scale(t(mydata),
                center = TRUE, #for cols (genes)
                scale = TRUE))

train <- mydata[rownames(mydata)%in%degs,] #78 290
mypdata$observations <- ifelse(mypdata$tissues=="Cancer",1, 0)

#-----------------------------------------------------------------#
# 2. Fitting lasso penalized logistic regression
#-----------------------------------------------------------------#

x <- t(train)
y <- mypdata$observations 

# Find the best lambda using cross-validation
# alpha=1 means lasso regression
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso, xlab= "log(λ)")+title("A. Cross-validated deviance of LASSO fit", line = 3)
abline(v=log(cv.lasso$lambda.min), col="blue", lty=5.5 )

# Compute the final model using the best lambda
lasso <- glmnet(x, y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)


lbs_fun <- function(lra, ...) {
  
  fit <- lra$glmnet.fit
  
  L=which(fit$lambda==lra$lambda.min)
  
  ystart <- sort(fit$beta[abs(fit$beta[,L])!=0,L])
  labs <- names(ystart)
  r <- range(fit$beta[,100]) # max gap between biggest and smallest coefs at smallest lambda i.e., 100th lambda
  yfin <- seq(r[1], r[2], length=length(ystart))
  
  xstart<- log(lra$lambda.min)
  xfin <- xstart+1
  
  text(xfin+0.3, yfin, labels=labs)
  segments(xstart, ystart, xfin, yfin)
}
plot(cv.lasso$glmnet.fit, xvar="lambda", label = TRUE, xlim=c(-8,0), lwd=2, xlab= "log(λ)")+title("A. Model coefficient paths", line = 3)
abline(v=log(cv.lasso$lambda.min), col="blue", lty=5.5 )
lbs_fun(cv.lasso)
text(-5.5,3.5,"v = log(0.0095)")

# Using lambda.min as the best lambda, gives the following regression coefficients
coef(cv.lasso, cv.lasso$lambda.min)
select.ind <- which(coef(cv.lasso, cv.lasso$lambda.min)!=0)
select.ind <- select.ind[-1]-1
select.variables <- rownames(train)[select.ind]
select.variables
# write.csv(select.variables, "./02_table/17_DEGs_LASSO.csv")


# Evaluation
## https://mp.weixin.qq.com/s/LCu0gtLGfy--E8mYXJaKuw
lasso.prob <- predict(lasso, newx = x, s = c(cv.lasso$lambda.min, cv.lasso$lambda.1se))
re <- cbind(y, lasso.prob)

library("ROCR")
library("caret")
pred_min <- prediction(re[,2],re[,1])
auc_min <- performance(pred_min, "auc")@y.values[[1]]
perf_min <- performance(pred_min, "tpr", "fpr")
plot(perf_min, colorize = FALSE, col = "blue")

library("ggplot2")
tpr_min = performance(pred_min,"tpr")@y.values[[1]]
dat = data.frame(tpr_min = perf_min@y.values[[1]],
                 fpr_min = perf_min@x.values[[1]])

ggplot() + 
    geom_line(data = dat, aes(x = fpr_min, y = tpr_min),color = "blue") + 
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
    theme_bw()+
    ggtitle("ROC curve")+
    annotate("text",x = .75, y = .25,
             label = paste("AUC = ",round(auc_min,2)),color = "blue")+
    scale_x_continuous(name  = "False positive rate (fpr)")+
    scale_y_continuous(name = "True positive rate (tpr)")


## 畫圖
library("png")
png("./03_figure/LASSO.png",    # create PNG for the heat map        
    width = 12*800,        # 5 x 300 pixels
    height = 10*800,
    res = 800)

par(mfrow = c(2,2))
plot(cv.lasso, xlab= "log(λ)")+title("A. Cross-validated deviance of LASSO fit", line = 3)
abline(v=log(cv.lasso$lambda.min), col="blue", lty=5.5 )

plot(cv.lasso$glmnet.fit, xvar="lambda", label = TRUE, xlim=c(-8,0), lwd=2, xlab= "log(λ)")+title("B. Model coefficient paths", line = 3)
abline(v=log(cv.lasso$lambda.min), col="blue", lty=5.5 )
lbs_fun(cv.lasso)
text(-6,3.5,"v = log(0.0095)")

plot(perf_min, colorize = FALSE, col = "red") + title("C. ROC curve for GEO datasets", line = 3)
abline(0,1)
text(0.7,0.2, paste("AUC = ",round(auc_min,2)))

img <- readPNG("C:/Users/flora/Desktop/gynecological cancers/03_figure/image.png")
rimg <- as.raster(img) # raster multilayer object
r <- nrow(rimg) / ncol(rimg) # image ratio
plot(NA, xlim = c(0, 1), ylim = c(0, r), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp=1, bty="n") + title("D. PPI network of 20 feature genes", line = 3)
rasterImage(rimg, 0, 0, 1+2/r, r)

dev.off()
