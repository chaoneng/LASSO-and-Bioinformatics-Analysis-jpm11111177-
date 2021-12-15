##' 2021-04-25
##' Identification of DE genes between cancer and normal tissues

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

source("./_function/DEGanalysis.R")
#-----------------------------------------------------------------#
# 1. Inputting three data sets of cervical cancer:
#    GSE9750, GSE7803, GSE63514
#-----------------------------------------------------------------#

df <- read.csv("./02_table/exprSet_CC_after_combat.txt", header = TRUE, sep = "\t", row.names = 1)
pdata <- read.csv("./02_table/Clinical_information_CC_GSE9750_GSE7803_GSE63514.csv", header = TRUE)

DEG_anlysis(df, pdata, "CC")

#-----------------------------------------------------------------#
# 2. Inputting three data sets of endometrial cancer:
#    GSE17025, GSE115810, GSE36389
#-----------------------------------------------------------------#

df <- read.csv("./02_table/exprSet_EC_after_combat.txt", header = TRUE, sep = "\t", row.names = 1)
pdata <- read.csv("./02_table/Clinical_information_EC_GSE17025_GSE115810_GSE36389.csv", header = TRUE)

DEG_anlysis(df, pdata, "EC")
