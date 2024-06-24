## Cluster AUCell
## Zhiyuan 
## 16 March 2023
## last modified 16 March 2023

library(data.table)
library(SeuratDisk)
library(SeuratData)
library(dplyr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(scales)
library(ggplot2)
stage_cols <- RColorBrewer::brewer.pal("Spectral", n = 11)
genotype_cols <- RColorBrewer::brewer.pal("Set2", n = 8)[c(1,2,3,8)]
sample_cols <- RColorBrewer::brewer.pal("Accent", n = 8)
doublet_cols <- RColorBrewer::brewer.pal("Paired", n = 8)[c(6,8,2)]
col.ls <- c(ArchR::ArchRPalettes$stallion, ArchR::ArchRPalettes$bear) # archR color palettes
names(col.ls) <- NULL
work_dir = 'GRN_scenicplus/ncall_newmotif/output/'


seu <- LoadH5Seurat("velocity/data/seurat/seu_RNAsoupx_NC.h5Seurat")
DefaultAssay(seu) <- "RNA"
seu2 <- seu[,seu$cell_type %in% c("NPB_nohox","NPB_nohox_cycling","dNC_nohox","dNC_nohox_cycling","mNC_nohox","mNC_head_mesenchymal","mNC_arch1","Pigment_gch2_high","Pigment_sox6_high")]
dim(seu2)

# filtered eRegulon
regulon_filtered <- fread(paste0(work_dir, "scenicplus/cytoscape/eRegulon_metadata_filtered.csv"))
df_regulon <- unique(regulon_filtered[,c(2,3,4,5)])
dim(df_regulon)
df_regulon$weird_names_rb <- gsub(pattern = "[+]",replacement = ".", df_regulon$Region_signature_name)
df_regulon$weird_names_rb <- gsub(pattern = "-",replacement = ".", df_regulon$weird_names_rb)
df_regulon$weird_names_rb <- gsub(pattern = "[(]",replacement = ".", df_regulon$weird_names_rb)
df_regulon$weird_names_rb <- gsub(pattern = "[)]",replacement = ".", df_regulon$weird_names_rb)
head(df_regulon$weird_names_rb)

df_regulon$weird_names_gb <- gsub(pattern = "[+]",replacement = ".", df_regulon$Gene_signature_name)
df_regulon$weird_names_gb <- gsub(pattern = "-",replacement = ".", df_regulon$weird_names_gb)
df_regulon$weird_names_gb <- gsub(pattern = "[(]",replacement = ".", df_regulon$weird_names_gb)
df_regulon$weird_names_gb <- gsub(pattern = "[)]",replacement = ".", df_regulon$weird_names_gb)
head(df_regulon$weird_names_gb)


# Read the AUC matrix and add it to the Seurat object
auc <- fread(paste0(work_dir, 'scenicplus/auc/eRegulon_AUC_region_based.tsv.gz'))
auc$barcode <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[1])
auc$sample <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[2])
auc$cell_id <- paste0(auc$sample, "_", auc$barcode)
auc <- data.frame(auc)
auc_mat <- auc[,2:936]
rownames(auc_mat) <- auc$cell_id
auc_mat <- as.matrix(auc_mat)
auc_mat <- t(auc_mat)
auc_mat <- auc_mat[,match(colnames(seu), colnames(auc_mat))]

auc_mat <- auc_mat[rownames(auc_mat) %in% df_regulon$weird_names_rb,]
rownames(auc_mat) <- df_regulon$Region_signature_name[match(rownames(auc_mat), df_regulon$weird_names_rb)]
auc_mat <- t(auc_mat)
write.csv(auc_mat, file = paste0(work_dir,"analysis/eRegulon_auc_rb_filtered.csv"), row.names = FALSE)
# seu[["eRegulon_AUC_region_based"]] <- CreateAssayObject(auc_mat)

auc_filtered <- auc_mat[rownames(auc_mat) %in% colnames(seu2),]
write.csv(auc_mat[rownames(auc_mat) %in% colnames(seu2),], file = paste0(work_dir,"analysis/eRegulon_auc_rb_filtered_wtnohox.csv"), row.names = FALSE)

## Read gene-based auc
auc <- fread(paste0(work_dir, 'scenicplus/auc/eRegulon_AUC_gene_based.tsv.gz'))
auc$barcode <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[1])
auc$sample <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[2])
auc$cell_id <- paste0(auc$sample, "_", auc$barcode)
auc <- data.frame(auc)
auc_mat <- auc[,2:936]
rownames(auc_mat) <- auc$cell_id
auc_mat <- as.matrix(auc_mat)
auc_mat <- t(auc_mat)
auc_mat <- auc_mat[,match(colnames(seu), colnames(auc_mat))]

auc_mat <- auc_mat[rownames(auc_mat) %in% df_regulon$weird_names_gb,]
rownames(auc_mat) <- df_regulon$Gene_signature_name[match(rownames(auc_mat), df_regulon$weird_names_gb)]
auc_mat <- t(auc_mat)
write.csv(auc_mat[rownames(auc_mat) %in% colnames(seu2),], file = paste0(work_dir,"analysis/eRegulon_auc_gb_filtered_wtnohox.csv"), row.names = FALSE)

