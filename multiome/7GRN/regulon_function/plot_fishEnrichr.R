# Plot FishEnrichr results
# Zhiyuan Hu
# 16 Jan 2024
# last modified 16 Jan 2024

# Load necessary library
library(ggplot2)

setwd("regulon_analysis/data/fishenrichr")
# Read the data

plist <- list()
for(itor in 0:9){
  df <- read.table(paste0("cluster",itor,"/GO_Biological_Process_2018_table.txt"), header = TRUE, sep = "\t")
  
  # Convert Adjusted P-value to a factor
  df$Adjusted_P_value <- cut(df$Adjusted.P.value, breaks = c(0, 0.001, 0.01, 0.05, 1), labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"))
  # Order the dataframe by Combined.Score and get the top 20
  if(nrow(df) > 20){
    df_top20 <- df[order(-df$Combined.Score), ][1:20, ]
  } else {
    df_top20 <- df
  }
  
  # Create the bar plot
  plist[[(itor+1)]] <- ggplot(df_top20, aes(y = Combined.Score, x = reorder(Term, Combined.Score), fill = Adjusted_P_value)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d() +  # Use a diverging color scale (you can change this)
    labs(y = "Combined Score", x = "Term", fill = "Adjusted P-value") +
    theme_classic() +
    coord_flip()  # Flips the axes
  
}
cowplot::plot_grid(plotlist = plist, ncol = 2)
ggsave('../../figures/fishenrichr/goer_barplots.pdf', width = 20, height = 30)
