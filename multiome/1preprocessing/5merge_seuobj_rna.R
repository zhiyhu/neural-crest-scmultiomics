# Gather all seurat object and make a tidy seurat obj for clustering
# Zhiyuan Hu
# 26 Dec 2022
# last modified 26 Dec 2022

library(Seurat)
singlet <- readRDS("multiome/analysis_newref/preprocessing/rds/rna_atac_singlet/metadata.rds")
head(singlet)

seu_list <- list()
for(i in 1:8){
  seu_list[[i]] <- readRDS(paste0("multiome/analysis_newref/preprocessing/rds/intermediate_seuobj/post_initialqc_scmo_s",i,"_RNA.rds"))
}

seu <- merge(seu_list[[1]], 
             y = c(seu_list[[2]], seu_list[[3]], seu_list[[4]], seu_list[[5]], seu_list[[6]], seu_list[[7]], seu_list[[8]]), 
             project = "foxd3_multiome")
rm(seu_list)

seu <- seu[,colnames(seu) %in% rownames(singlet)]
singlet <- singlet[colnames(seu),]
seu@meta.data <- cbind(seu@meta.data, singlet[,2:17])

saveRDS(seu, "multiome/analysis_newref/preprocessing/rds/rna_atac_singlet/seuobj_rna64340.rds")
