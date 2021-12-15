## 2021-04-25
## Creating a heat map of the DE genes

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

#-----------------------------------------------------------------#
# Data input
#-----------------------------------------------------------------#

library("dplyr")

filenames <- list.files(path = "./02_table/", pattern = "_after_combat.txt", full.names=TRUE)
All <- lapply(filenames, function(i){
  read.csv(i, header = TRUE, sep = "\t", row.names = 1)
})
names(All) <- c("CC", "CC_EC", "EC")

DEG_CC <- read.csv("./02_table/DEG_CC_Limma.csv", header = TRUE, sep = ",", row.names = 1)
DEG_EC <- read.csv("./02_table/DEG_EC_Limma.csv", header = TRUE, sep = ",", row.names = 1)

x <- list(
  DEG_CC %>% filter(abs(logFC)>1&adj.P.Val<0.05) %>% rownames() %>% unlist(),
  DEG_EC %>% filter(abs(logFC)>1&adj.P.Val<0.05) %>% rownames() %>% unlist())
common <- intersect(x[[1]], x[[2]])

#-----------------------------------------------------------------#
# 2. Heat map
#-----------------------------------------------------------------#

source("./_function/plotHeamap.R")

CC <- All[["CC"]]
mypdata <- read.csv("./02_table/Clinical_information_CC_GSE9750_GSE7803_GSE63514.csv", header = TRUE, sep = ",", row.names = 1)
heamap_plot(CC, mypdata, common, "CC")
heamap_plot(CC, mypdata, degs, "CC")

EC <- All[["EC"]]
colnames(EC) <- sapply(strsplit(colnames(EC), "_"), "[", 1)
mypdata2 <- read.csv("./02_table/Clinical_information_EC_GSE17025_GSE115810_GSE36389.csv", header = TRUE, sep = ",", row.names = 1)
heamap_plot(EC, mypdata2, common, "EC")
heamap_plot(EC, mypdata2, degs, "EC")

CC_EC <- All[["CC_EC"]]
mypdata3 <- read.csv("./02_table/Clinical_information_CC_EC.csv", header = TRUE, sep = ",", row.names = 1)
degs <- read.csv("./02_table/16_DEGs_LASSO.csv", row.names = 1)
degs <- degs$x 

heamap_plot(mydata, mypdata3, degs, "CC_EC")
