##' 2021-04-26
##' Combat()
##' Ref: https://blog.csdn.net/yayiling/article/details/113664434

rm(list = ls())
Sys.setenv(LANGUAGE="en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

library("sva")

cc <- read.table("./02_table/exprSet_CC_before_combat.txt", header = TRUE, sep = "\t", row.names = 1)
ec <- read.table("./02_table/exprSet_EC_before_combat.txt", header = TRUE, sep = "\t", row.names = 1)

pdata <- read.csv("./02_table/Clinical_information_CC_EC.csv", header = TRUE)

mrna_names <- intersect(rownames(cc),rownames(ec))
exprSet <- cbind(cc[mrna_names,], ec[mrna_names,])
write.table(exprSet, "./02_table/exprSet_CC_EC_before_combat.txt", sep = "\t", col.names = NA)

batch <- paste0("batch", rep(c(1,2), c(140,150)))
tissue <- pdata$tissues
table(batch, tissue)
mod <- model.matrix(~tissue)

#---- Using Combat()

exprSet_batch <- ComBat(dat = exprSet, batch = batch, mod = mod)
write.table(exprSet_batch, "./02_table/exprSet_CC_EC_after_combat.txt", sep = "\t", col.names = NA)
