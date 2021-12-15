
require("ggplot2")
require("gplots")
require("pheatmap") 
require("ComplexHeatmap")
require("circlize")
require("RColorBrewer")
library("IMvigor210CoreBiologies")

heamap_plot <- function(expr, pdata, commonGenes, type){
  
  if(type=="CC"){
    fig_title <- paste0("Heatmap of the ", length(commonGenes), " DEGs in cervical cancer")
    png_path <- paste0("./03_figure/", fig_title, ".png")
  }else if(type=="EC"){
    fig_title <- paste0("Heatmap of the ", length(commonGenes), " DEGs in endometrial cancer")
    png_path <- paste0("./03_figure/", fig_title, ".png")
  }else{
    fig_title <- paste0("Heatmap of the ", length(commonGenes), " DEGs in cervical and endometrial cancers")
    png_path <- paste0("./03_figure/", fig_title, ".png")
    
  }
  heat <- t(scale(t(expr)))
  heat <- heat[rownames(heat)%in%commonGenes,]

  ann_colors <- list("Disease state" = c(Cancer="Red", Normal="#027CD9"))
  
  ha <- HeatmapAnnotation("Disease state"=pdata$tissues,
                          annotation_legend_param=list(labels_gp = gpar(fontsize = 9),
                                                       title_gp = gpar(fontsize = 9, fontface = "bold"),
                                                       ncol=1),
                          gap=unit(c(1, rep(0, 4), 1, rep(0, 5)), "mm"),
                          col=list("Disease state"=ann_colors$`Disease state`),
                          show_annotation_name = TRUE,
                          show_legend = FALSE,
                          annotation_name_side="right",
                          annotation_name_gp = gpar(fontsize = 10, fontface = "bold"))
  
  heat_colors2 <- colorRamp2(c(-2.5, 0, 2.5), c("green", "black", "red"))
  
  htmap <- Heatmap(limitRange(heat), ##對基因表現量進行限制
                   name="Expression",
                   top_annotation = ha,#頂部注釋
                   cluster_rows = TRUE,#對rows進行聚類
                   cluster_columns = TRUE, #對cols進行聚類
                   col=heat_colors2,#基因表現量顏色
                   color_space = "RGB",
                   row_order=NULL,
                   column_order=NULL,
                   show_column_names = FALSE,
                   #show_row_names = FALSE,
                   row_names_gp = gpar(fontsize = 8),
                   gap = unit(1, "mm"),
                   column_title = fig_title,
                   column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                   width=unit(10, "cm"),
                   show_heatmap_legend = FALSE,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 9), 
                                             title_gp = gpar(fontsize = 9, fontface = "bold"))
  )
  ht_list <- htmap
  
  lgd1 <- Legend(labels = names(ann_colors$`Disease state`), title = "Disease state",
                 legend_gp = gpar(fill = ann_colors$`Disease state`),
                 nr = 2)
  # The discrete legend for continuous color mapping:
  at = seq(-4, 4, by = 2)
  lgd2 <- Legend(at = at, title = "Expression", 
                 legend_gp = gpar(fill = heat_colors2(at)))
  
  draw(ht_list, 
       padding = unit(c(2,0,7,0), "mm"), 
       heatmap_legend_list = list(lgd1, lgd2))
  
  # Saving a png
  png(png_path,           
      width = 7*900,        # 5 x 300 pixels
      height = 8*900,
      res = 900)             # 300 pixels per inch 
  
  draw(ht_list, 
       padding = unit(c(2,0,7,0), "mm"), 
       heatmap_legend_list = list(lgd1, lgd2))
  
  dev.off()
}