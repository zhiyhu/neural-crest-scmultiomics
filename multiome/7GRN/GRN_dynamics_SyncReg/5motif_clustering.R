#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cluster the motifs based on TOMTOM similarity
# Authors: Zhiyuan Hu, ChatCPT
# Date: 25 Apr 2023
# Last modified 25 Apr 2023
# refer to: https://github.com/bernardo-de-almeida/motif-clustering/blob/main/Motif_clustering.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the below script first
# latent_time_grn/scripts/tomtom.sh
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(fields)
setwd("latent_time_grn/figures/motif_clustering/")

# read the motif2TF metadata-----
metadata <- read.delim("GRN_scenicplus/data/motif_with_newm/dr11_motif2tf_curated_filtered20230122.tbl", as.is = TRUE)
metadata <- metadata[,-1]
colnames(metadata)[1] <- "motif_id"
head(metadata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Hierarchically cluster motifs by similarity ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# based on Jeff Vierstra https://raw.githubusercontent.com/jvierstra/motif-clustering/master/hierarchical.py

## Load similarity data ----
sim_file <- "GRN_scenicplus/data/motif_with_newm/TOMTOM/tomtom.all.treated.txt"
sim <- read.delim(sim_file)

## Transform similarity data into a square matrix -----
# E-value: multiplying the p-value by the total number of target motifs in all the target databases
simsq <- data.table::dcast(sim, Query_ID~Target_ID, value.var = "E.value")
rownames(simsq) <- simsq$Query_ID
simsq <- simsq[,-1]
simsq[is.na(simsq)] <- 100
simsq[1:20,1:20]

# Calculate the Pearson correlation matrix from the similarity matrix
mat = -log10(simsq)
mat[mat == Inf]=10
tmp_cor <- cor(mat, method="pearson")

# Perform hierarchical clustering using Pearson correlation as the distance metric and complete linkage
Z = hclust(as.dist(1-tmp_cor), method = 'complete')

## Test clustering at a range of tree heights (0.5-1) -----
thresholds=seq(0.5, 1.0, 0.1)

pdf("Hierarchical_clusters_diff_thresholds.pdf", width = 20, height = 5)
for(thresh in thresholds){
  
  cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
  df = data.frame(Motifs=names(mat),
                  Cluster=cl[match(names(mat), names(cl))])
  write.table(df, paste0("clusters.", thresh,".txt"), sep="\t", row.names = F, quote=F)
  
  plot(Z, labels=FALSE, main=paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))
  abline(h=thresh, col="red")
  
  print(paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))

}
dev.off()

## Choose final clusters cutting the dendrogram at height 0.8 -----
thresh=0.8
cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
df = data.frame(Motifs=names(mat),
                Cluster=cl[match(names(mat), names(cl))],
                Order_dendogram=match(names(mat), Z$labels[Z$order]))
# df <- merge(df, TF_clusters_PWMs$metadata[,c(1,13,10)], by=1)
df <- df[order(df$Order_dendogram),]
length(unique(df$Cluster))
sort(table(df$Cluster))
save(Z, mat, file = paste0("All_motifs_data_and_hclust_objects.Rdata"))
saveRDS(df, paste0("All_motifs_final_clusters_thresh", thresh, ".rds"))

## plot hierarchical clustering heatmap of motifs clustered by simililarity ----
##  and  clusters identified cutting the dendrogram at height 0.8
# Notice that image interprets the z matrix as a table of f(x[i], y[j]) values, 
# so that the x axis corresponds to row number and the y axis to column number,
# with column 1 at the bottom, i.e. a 90 degree counter-clockwise rotation of 
# the conventional printed layout of a matrix.
# that's why I need to reverse the order of the columns
# top-right should represent cluster1. - and so on
png(paste0("All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png"), 
    type="cairo", width = 2200, height = 2000, res = 300)
out <- tmp_cor[Z$order,rev(Z$order)]
# Create a custom color gradient from blue to green to red
color_gradient <- colorRampPalette(c("white", "blue"))(100)
image.plot(out, col=color_gradient, #col=c("white", "black"), # colorRampPalette(c("grey100", "grey0"))(100)
      legend.lab="Value", legend.width=2, las=1, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
dev.off()

## Highlight specific clusters on the hierarchical clustering heatmap ----
for(c in 1:length(unique(df$Cluster))){
  png(paste0("Highlight_cluster_",c,".png"), type="cairo", width = 2000, height = 2000, res = 50)
  highli <- names(cl)[which(cl==c)]
  image.plot(out, col=color_gradient,
        las=1, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
  axis(1, at=seq(0,1,length.out = nrow(out))[rownames(out) %in% highli], 
       labels=rep("",length(highli)), col = "red")
  axis(2, at=seq(0,1,length.out = nrow(out))[colnames(out) %in% highli], 
       labels=rep("",length(highli)), col = "red")
  dev.off()
  print(c)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Integrate clustering results into the metadata ----
metadata$motif_clusters <- NA
clusters_df <- read.delim("clusters.0.8.txt")
metadata$motif_clusters <- clusters_df$Cluster[match(metadata$motif_id, clusters_df$Motifs)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subset the similarity matrix by the Regulons included in our analysis ----
## Plot heatmap and dendrogram
auc_df <- read.csv("regulon_analysis/figures/process_regulon_data/auc_gb_clean.csv", as.is = TRUE, header = FALSE)
regulons_gb <- unlist(auc_df[1,])
auc_df <- read.csv("regulon_analysis/figures/process_regulon_data/auc_rb_clean.csv", as.is = TRUE, header = FALSE)
regulons_rb <- unlist(auc_df[1,])

TFs <- unique(c(sapply(regulons_gb, function(x) unlist(strsplit(x, "_"))[1]),
                     sapply(regulons_rb, function(x) unlist(strsplit(x, "_"))[1])))
head(TFs)

metadata_subset <- metadata[metadata$gene_name %in% TFs, ]
tmp_cor_subset <- tmp_cor[metadata_subset$motif_id,metadata_subset$motif_id ]
rownames(tmp_cor_subset) <- colnames(tmp_cor_subset) <- metadata_subset$gene_name
rownames(tmp_cor_subset)[grepl("glam", metadata_subset$motif_id)] <- colnames(tmp_cor_subset)[grepl("glam", metadata_subset$motif_id)] <- "foxd3_glam2"
Z = hclust(as.dist(1-tmp_cor_subset), method = 'complete')

## Plot the dendrogram with leaf names----
pdf("dendrogram_of_motifs_TFnames.pdf", width = 12, height = 8)
plot(Z, labels = rownames(tmp_cor_subset) , cex = 0.6, xlab = "Motifs", ylab = "Height", main = "Dendrogram of Motifs")
dev.off()

# Write the tmp_cor data
write.table(tmp_cor_subset, file = "similarity_matrix_subset.csv", sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write session information to the file -----
session_info_file <- "latent_time_grn/sessionInfo/motif_clustering_sessionInfo.txt"
session_info_output <- capture.output(sessionInfo())
writeLines(session_info_output, session_info_file)
