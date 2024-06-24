## preproces regulon + latent time data
## Zhiyuan Hu
## 22 APR 2023
## last modified 17 Jan 2024

library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(gridExtra)

setwd("latent_time_grn")
work_dir = 'GRN_scenicplus/ncall_2023oct_ccb/output/'

mv_time <- read.csv("data/wtnohox/multivelo_obs.csv")
scv_time <- read.csv("data/wtnohox/scvelo_obs.csv")
head(scv_time)

# trajectories
traj1 <- c( "mNC_head_mesenchymal", "mNC_nohox","dNC_nohox_cycling", "dNC_nohox", "NPB_nohox_cycling", "NPB_nohox")
traj2 <- c( "Pigment_gch2_high", "Pigment_sox6_high", "mNC_nohox", "dNC_nohox_cycling", "dNC_nohox", "NPB_nohox_cycling", "NPB_nohox")
traj3 <- c( "mNC_arch1", "mNC_nohox","dNC_nohox_cycling" ,  "dNC_nohox",   "NPB_nohox_cycling", "NPB_nohox" )

# keep all trajectories
idx <- mv_time$X[mv_time$cell_type %in% c(traj1, traj2, traj3)]

# regulon set clustering results
regulon_sets <- read.csv("regulon_analysis/figures/cluster_regulons/clustering_results_gb.csv")

##~~~~~~~~~~~~~~~~~~~##
## Region based
##~~~~~~~~~~~~~~~~~~~##
auc <- fread(paste0(work_dir, 'scenicplus/auc/eRegulon_AUC_gene_based.tsv.gz'))
auc[1:5,1:5]
auc_colnames <- colnames(auc)
auc$barcode <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[1])
auc$sample <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[2])
auc$cell_id <- paste0(auc$sample, "_", auc$barcode)
# auc <- data.frame(auc)
auc[1:5,1:5]
auc_cells <- auc$cell_id
auc <- auc[, 2:(ncol(auc)-3)]
auc_regulons <- colnames(auc)
auc <- as.data.frame(auc)


# plot PDF scatter plots for each set
for(itor_set in unique(regulon_sets$Agglomerative)){
 set_n <-  regulon_sets$Regulon[regulon_sets$Agglomerative == itor_set]
 set_n 
  
  # Create the data frame
  plist <- list()
  for(i in 1:length(set_n)){
    df_plot <- data.frame(cell_id = idx,
                          latent_time = mv_time$latent_time[match(idx, mv_time$X)],
                          cell_type = mv_time$cell_type[match(idx, mv_time$X)],
                          AUC_scores = auc[match(idx, auc_cells), auc_regulons ==set_n[i]])
    head(df_plot)
    # Plot the data
    plist[[i]] <- ggplot(df_plot, aes(x = latent_time, y = AUC_scores, col = cell_type)) +
      geom_point() +
      labs(title =set_n[i]) +
      theme_light()
  }
  
  # save the plot
  pdf(paste0("latent_time_grn/figures/process_regulon_data/scatterplots_gb_",itor_set,".pdf"), onefile = TRUE, width = 5, height = 3)
  for (i in seq(length(plist))) {
    print(plist[[i]])  
  }
  dev.off()
}

## plot the pngs for each of the regulon
for(i in 1:length(regulon_sets$Regulon)){

  itor_regulon <- regulon_sets$Regulon[i]
  df_plot <- data.frame(cell_id = idx,
                        latent_time = mv_time$latent_time[match(idx, mv_time$X)],
                        cell_type = mv_time$cell_type[match(idx, mv_time$X)],
                        AUC_scores = auc[match(idx, auc_cells), colnames(auc) == itor_regulon])

    p <- ggplot(df_plot, aes(x = latent_time, y = AUC_scores, col = cell_type)) +
      geom_point() +
      theme_void() + 
      theme(legend.position = "none")
  ggsave(plot = p, filename = paste0("latent_time_grn/figures/regulon_auc_gb_scatterplots/",itor_regulon,".png"), width = 6, height = 6)
}

# save the data into csv files
auc_df <- auc[match(idx, auc_cells), auc_regulons %in% regulon_sets$Regulon]
write.csv(x = auc_df, "latent_time_grn/figures/process_regulon_data/auc_gb_clean.csv", row.names = FALSE)

