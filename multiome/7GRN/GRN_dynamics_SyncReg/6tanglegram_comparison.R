# Script to Compare Hierarchical Clustering Results
# Author: Zhiyuan, ChatGPT
# Date: 2023-04-25
# Last modified: 2023-04-25
#
# This script performs hierarchical clustering on a similarity matrix and compares the resulting dendrograms
# with two other dendrograms from different similarity matrices using tanglegrams.
#
library(dendextend )
setwd("latent_time_grn/figures/")
# Function to perform hierarchical clustering and create a dendrogram
perform_clustering_dendrogram <- function(similarities) {
  # Calculate distance matrix
  dist_matrix <- as.dist(1 - similarities)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = "complete")
  
  # Create a dendrogram object
  dend <- as.dendrogram(hc)
  
  return(dend)
}

# Perform hierarchical clustering and create a dendrogram for the original similarity matrix
motif_sim <-  read.csv("motif_clustering/similarity_matrix_subset.csv", row.names = 1)

sim_mat1 <- read.csv("image_similarity_rb/similarity_matrix.csv", row.names = 1)
sim_mat2 <- read.csv("image_similarity_gb/similarity_matrix.csv", row.names = 1)
colnames(sim_mat1) <- rownames(sim_mat1) <- gsub("_\\(.*\\)", "", rownames(sim_mat1))
colnames(sim_mat2) <- rownames(sim_mat2) <- gsub("_\\(.*\\)", "", rownames(sim_mat2))

dendrb <- perform_clustering_dendrogram(sim_mat1)
dendgb <- perform_clustering_dendrogram(sim_mat2)

### Region based vs Motif based
tf_df <- data.frame(rb_regulon = rownames(sim_mat1),
                    tf = sapply(rownames(sim_mat1), function(x) unlist(strsplit(x, "_"))[1]))

idx <- match(tf_df$tf, rownames(motif_sim))
motif_sim_modified <- motif_sim[idx, idx]
rownames(motif_sim_modified) <- colnames(motif_sim_modified) <- tf_df$rb_regulon
motif_dend <- perform_clustering_dendrogram(motif_sim_modified)
labels(motif_dend)

pdf("compare_dend/tanglegram_comparison_motifBased_regionBased.pdf", width = 15, height = 20)
tanglegram(motif_dend, dendrb,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_branches = TRUE, # Color common branches
           columns_width = c(6,3,6),
           left_dendo_mar =  c(3, 3, 3, 6),
           right_dendo_mar = c(3, 6, 3, 3),
           main_left = "Motif-based",
           main_right = "Region-based")
dev.off()

### Gene based vs Motif based
tf_df <- data.frame(gb_regulon = rownames(sim_mat2),
                    tf = sapply(rownames(sim_mat2), function(x) unlist(strsplit(x, "_"))[1]))

idx <- match(tf_df$tf, rownames(motif_sim))
motif_sim_modified <- motif_sim[idx, idx]

rownames(motif_sim_modified) <- colnames(motif_sim_modified) <- tf_df$gb_regulon
motif_dend <- perform_clustering_dendrogram(motif_sim_modified)
labels(motif_dend)

pdf("compare_dend/tanglegram_comparison_motifBased_geneBased.pdf", width = 15, height = 20)
tanglegram(motif_dend, dendgb,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_branches = TRUE, # Color common branches
           columns_width = c(6,3,6),
           left_dendo_mar =  c(3, 3, 3, 6),
           right_dendo_mar = c(3, 6, 3, 3),
           main_left = "Motif-based",
           main_right = "Gene-based")
dev.off()
