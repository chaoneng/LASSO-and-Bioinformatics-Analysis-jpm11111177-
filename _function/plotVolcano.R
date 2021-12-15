
require("ggplot2")
volcano_plot <- function(df, FC, FDR, type){
  
  if(type=="CC"){
    sub_title <- "cervical cancer vs. normal tissue"
  }else{
    sub_title <- "endometrial cancer vs. normal tissue"
  }
  
  DEG <- df %>% 
    mutate(labels = ifelse(df$logFC>FC&df$adj.P.Val<FDR, "up regulated", 
                          ifelse(df$logFC<(-FC)&df$adj.P.Val<FDR, "down regulated", "non-significant")))
  
  p1 <- ggplot(DEG, aes(x = logFC, y = -log10(adj.P.Val), colour = labels)) +
    labs(x = expression(log[2]("Fold Change")), y = expression(-log[10]("adj.P.Val")),
         title = "Volcano plot for differentially expressed genes",
         subtitle = sub_title) +
    theme_bw()+
    xlim(-2.0,2.0) +
    ylim(0.0,25.0)+
    geom_point(size = 1.5) +
    scale_color_manual(values = c("up regulated" = "red", "non-significant" = "black", "down regulated" = "green")) +
    geom_vline(xintercept = c(-log2(2),log2(2)), lty=2, size=I(0.5), colour = "blue") +
    geom_hline(yintercept = -log10(0.1), lty=2, size=I(0.5), colour = "blue")
  
}
