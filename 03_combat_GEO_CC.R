##' 2021-04-25
##' Combat()
##' Ref: https://blog.csdn.net/yayiling/article/details/113664434

rm(list = ls())
Sys.setenv(LANGUAGE="en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

library("sva")

listSAmples <- c("GSE9750", "GSE7803", "GSE63514",
                 "GSE17025", "GSE115810", "GSE36389")
pdata <- read.csv("./02_table/Clinical_information_CC_GSE9750_GSE7803_GSE63514.csv", header = TRUE)

#---- Loading required data sets

filenames <- NULL
for(i in listSAmples){
  filename <- paste0("./01_data/", i, "_rma_symbol.txt")
  filenames <- c(filenames, filename)
}
allSamples <- lapply(filenames, function(i){
  read.table(i, header = TRUE, sep = "\t", row.names = 1)
})
names(allSamples) <- listSAmples

# Cervical cancer: "GSE9750", "GSE7803", "GSE63514"
GSE9750 <- allSamples[["GSE9750"]]
GSE7803 <- allSamples[["GSE7803"]]
GSE63514 <- allSamples[["GSE63514"]]

#---- Intersection of three data sets from GEO

mrna_names <- intersect(rownames(GSE9750),rownames(GSE7803))
mrna_intersect <- intersect(mrna_names, rownames(GSE63514))
exprSet <- cbind(GSE9750[mrna_intersect,], GSE7803[mrna_intersect,], GSE63514[mrna_intersect,])
write.table(exprSet, "./02_table/exprSet_CC_before_combat.txt", sep = "\t", col.names = NA)

batch <- paste0("batch", rep(c(1,2,3), c(57,31,52)))
tissue <- pdata$tissues
table(batch, tissue)
mod <- model.matrix(~tissue)

#---- Using Combat()

exprSet_batch <- ComBat(dat = exprSet, batch = batch, mod = mod)
write.table(exprSet_batch, "./02_table/exprSet_CC_after_combat.txt", sep = "\t", col.names = NA)


library("gplots")
heatmap.2(as.matrix(exprSet[1:500,]),cexrow=0.8,cexcol=1.0)
heatmap.2(as.matrix(exprSet_batch[1:500,]),cexrow=0.8,cexcol=1.0)


