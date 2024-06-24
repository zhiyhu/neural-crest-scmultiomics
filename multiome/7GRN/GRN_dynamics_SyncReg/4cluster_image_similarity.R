# Image Similarity Heatmap using pheatmap
# Author: ChatGPT, Zhiyuan
# Date: 2023-04-20
# Last modified: 2023-04-25

# Load required libraries
library(readr)
library(pheatmap)
library(RColorBrewer)
library(dendextend)

setwd("latent_time_grn/figures/")

# Read and prepare data
prepare_data <- function(file_path) {
  data <- read_csv(file_path)
  data_mat <- data.frame(data, row.names = 1)
  colnames(data_mat) <- rownames(data_mat)
  return(data_mat)
}

# Plot similarity heatmap
plot_heatmap <- function(similarities_mat, output_file) {
  pheatmap(similarities_mat,
           color = colorRampPalette(c("darkblue", "white", "darkred"))(256),
           fontsize = 5,
           cellwidth = 4, cellheight = 4,
           border_color = NA,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           filename = output_file)
}

# Perform hierarchical clustering and create a dendrogram
perform_clustering_dendrogram <- function(similarities, method = "complete") {
  dist_matrix <- as.dist(1 - as.matrix(similarities))
  hc <- hclust(dist_matrix, method = method)
  dend <- as.dendrogram(hc)
  return(dend)
}

# Plot horizontal dendrogram
plot_dendrogram <- function(dend, output_file) {
  pdf(output_file, width = 8, height = 30)
  par(mar = c(3, 3, 3, 10) + 0.1) # Increase left margin to make room for the labels
  plot(dend, horiz = TRUE, main = "Hierarchical Clustering Dendrogram", nodePar = list(cex = 0.6))
  dev.off()
}

# Main script
setwd("latent_time_grn/figures/")

sim_mat1 <- prepare_data("image_similarity_rb/similarities_with_labels.csv")
sim_mat2 <- prepare_data("image_similarity_gb/similarities_with_labels.csv")

sim_mat3 <- prepare_data("image_similarity_colHist/similarities_with_labels.csv")

# Save similarity matrices as CSV files
write.csv(sim_mat1, "image_similarity_rb/similarity_matrix.csv", row.names = TRUE)
write.csv(sim_mat2, "image_similarity_gb/similarity_matrix.csv", row.names = TRUE)

write.csv(sim_mat3, "image_similarity_colHist/similarity_matrix.csv", row.names = TRUE)


plot_heatmap(sim_mat1, "image_similarity_rb/similarity_heatmap_pheatmap.pdf")
plot_heatmap(sim_mat2, "image_similarity_gb/similarity_heatmap_pheatmap.pdf")

plot_heatmap(sim_mat3, "image_similarity_colHist/similarity_heatmap_pheatmap.pdf")


# Perform hierarchical clustering and create dendrograms for both similarity matrices
colnames(sim_mat1) <- rownames(sim_mat1) <- gsub("_\\(.*\\)", "", rownames(sim_mat1))
colnames(sim_mat2) <- rownames(sim_mat2) <- gsub("_\\(.*\\)", "", rownames(sim_mat2))

colnames(sim_mat3) <- rownames(sim_mat3) <- gsub("_\\(.*\\)", "", rownames(sim_mat3))


dend1 <- perform_clustering_dendrogram(sim_mat1)
dend2 <- perform_clustering_dendrogram(sim_mat2)

dend3 <- perform_clustering_dendrogram(sim_mat3)


plot_dendrogram(dend1, "image_similarity_rb/horizontal_dendrogram.pdf")
plot_dendrogram(dend2, "image_similarity_gb/horizontal_dendrogram.pdf")

plot_dendrogram(dend3, "image_similarity_colHist/horizontal_dendrogram.pdf")


# Compare dendrograms using tanglegram
pdf("compare_dend/tanglegram_comparison_imgsim.pdf", width = 15, height = 20)
tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_branches = TRUE, # Color common branches 
           columns_width = c(6,3,6),
           left_dendo_mar =  c(3, 3, 3, 6),
           right_dendo_mar = c(3, 6, 3, 3),
           main_left = "Region-based",
           main_right = "Gene-based")
dev.off()
