##' 2021-04-27
##' Gene function analysis of 43 DE genes
##' https://www.jianshu.com/p/a55244438bd5

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
setwd("C:/Users/flora/Desktop/gynecological cancers") 

#--- Preparing required packages ------
library("clusterProfiler")
library("org.Hs.eg.db")
library("cowplot")
library("ggplot2")
library("stringr")

#--- Loading required data ------
degs <- read.csv("./02_table/78_common_DE_genes.csv", row.names = 1)
degs <- degs$x
g <- bitr(degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#--- Analyze and visualize functional profiles (GO and KEGG) ------
##' Gene ontology
go.all <- enrichGO(gene = g$ENTREZID,
                   OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                   #keytype = 'ENSEMBL',
                   ont = "ALL", #也可以是 CC  BP  MF中的一种
                   pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                   pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                   qvalueCutoff = 0.1,
                   readable = TRUE)
go.bar <- barplot(go.all, showCategory = 10, font.size = 10, title = "A. GO enrichment", split = "ONTOLOGY", color = "pvalue") + 
  facet_grid(ONTOLOGY ~ ., scales = "free")+
  scale_y_discrete(labels=function(y) stringr::str_wrap(y, width = 46))

# setReadable(go.all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

##' KEGG
enrich.kegg <- enrichKEGG(gene = g$ENTREZID,
                          organism ="hsa",
                          keyType = "kegg",
                          minGSSize = 2,
                          maxGSSize = 77,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          use_internal_data =FALSE)
kegg.bar <- barplot(enrich.kegg, showCategory=15, font.size = 10, title = "B. KEGG enrichment", color = "pvalue") + 
  xlim(NA,4) +
  scale_y_discrete(labels=function(y) stringr::str_wrap(y,width=34))

#--- Saving a plot ------
png("./03_figure/GO_KEGG_clusterProfiler.png",    # create PNG for the heat map        
    width = 16*800,        # 5 x 300 pixels
    height = 8*800,
    res = 800)

plot_grid(go.bar, kegg.bar, ncol=2, rel_widths = c(5.5,4))

dev.off()
