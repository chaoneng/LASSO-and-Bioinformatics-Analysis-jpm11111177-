##' 2021-04-24
##' Integrating samples information

rm(list = ls())
Sys.setenv(LANGUAGE="en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

listSAmples <- c("GSE9750", "GSE7803", "GSE63514",
                 "GSE17025", "GSE115810", "GSE36389")
filenames <- NULL
for(i in listSAmples){
  filename <- paste0("./01_data/", i, "_clinical.csv")
  filenames <- c(filenames, filename)
}
allSamples <- lapply(filenames, function(i){
  read.csv(i, header = TRUE, sep = ",", row.names = 1)
})
names(allSamples) <- listSAmples

index <- c("geo_accession", "tissues")

##' Cervical cancer: "GSE9750", "GSE7803", "GSE63514"

GSE9750 <- allSamples[["GSE9750"]]
GSE9750$tissues <- rep(c("Normal", "Cancer"), c(24, 33))
GSE9750 <- GSE9750[,index]

GSE7803 <- allSamples[["GSE7803"]]
GSE7803$tissues <- rep(c("Normal", "Cancer"), c(10, 21))
GSE7803 <- GSE7803[,index]

GSE63514 <- allSamples[["GSE63514"]]
GSE63514$tissues <- rep(c("Normal", "Cancer"), c(24, 28))
GSE63514 <- GSE63514[,index]

CC_groups <- rbind(GSE9750, GSE7803)
CC_groups <- rbind(CC_groups, GSE63514)
CC_groups$GEO <- rep(c("GSE9750", "GSE7803", "GSE63514"), c(57,31,52)) 
# write.csv(CC_groups, "./02_table/Clinical_information_CC_GSE9750_GSE7803_GSE63514.csv", row.names = FALSE)

##' Endometrial cancer: "GSE17025", "GSE115810", "GSE36389"

GSE17025 <- allSamples[["GSE17025"]]
GSE17025$tissues <- rep(c("Cancer", "Normal"), c(91, 12))
GSE17025 <- GSE17025[,index]

GSE115810 <- allSamples[["GSE115810"]]
GSE115810$tissues <- rep(c("Cancer", "Normal", "Cancer"), c(2,3, 22))
GSE115810 <- GSE115810[,index]

GSE36389 <- allSamples[["GSE36389"]]
GSE36389$tissues <- ifelse(GSE36389$cell.type.ch1=="endometrium control", "Normal", "Cancer")
GSE36389 <- GSE36389[,index]

EC_group <- rbind(GSE17025, GSE115810)
EC_group <- rbind(EC_group, GSE36389)
EC_group$GEO <- rep(c("GSE17025", "GSE115810", "GSE36389"), c(103,27,20))
# write.csv(EC_group, "./02_table/Clinical_information_EC_GSE17025_GSE115810_GSE36389.csv", row.names = FALSE)
