# Integrate Multiome data all cells with Wagner 2018
# Zhiyuan Hu
# 15 Dec 2022
# last modified 27 Dec 2022
## refer to https://satijalab.org/seurat/articles/integration_large_datasets.html

## prepare multiome data


library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(loomR)
library(ggplot2)
library(cowplot)
library(patchwork)

figdir="/t1-data/project/tsslab/zhu/multiome/analysis_newref/integration/figures/allcells_wagner_cca/"
outputdir="/t1-data/project/tsslab/zhu/multiome/analysis_newref/integration/data/"
ndim=50
##---------------------------------##
## Preprocessing both datasets ------
##---------------------------------##
seu <- readRDS("/t1-data/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NTNC31037.rds")

seu <- CreateSeuratObject(counts = seu@assays$RNA@counts,
                          meta.data = seu@meta.data)
seu

# SaveH5Seurat(seu, filename = "/t1-data/project/tsslab/zhu/multiome/R/clustering/rds/subset_NTnNC/NTnNC_35241.h5Seurat")
# Convert("/t1-data/project/tsslab/zhu/multiome/R/clustering/rds/subset_NTnNC/NTnNC_35241.h5Seurat", dest = "h5ad")

## preprocessing Wagner 2018 data------
# https://satijalab.org/loomr/loomr_tutorial
refloom <- connect(filename = "/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.loom", mode = "r", skip.validate = TRUE)
# ref <- Convert(refloom, to = "seurat")

full.matrix <- refloom[["matrix"]][,]
full.matrix <- t(full.matrix)
dim(x = full.matrix)

# Access all gene names
gene.names <- refloom[["row_attrs/index"]][]
rownames(full.matrix) <- gene.names

# Access all cell names
cell.names <- refloom$col.attrs$index[]
colnames(full.matrix) <- cell.names
gc()

# metadata
ref_metadata <- read.csv("/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018_metadata.csv", row.names = 1)

# create seurat object
ref <- CreateSeuratObject(counts = full.matrix, meta.data = ref_metadata)
rm(full.matrix)
rm(ref_metadata)
refloom$close_all()

ref$ClusterName_short <- sapply(ref$ClusterName, function(x) unlist(strsplit(x, "hpf-"))[2])
# SaveH5Seurat(ref, filename = "/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.processed.h5Seurat")
# Convert("/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.processed.h5Seurat", dest = "h5ad")


## 16755 gene overlap
table(rownames(ref) %in% rownames(seu))
# FALSE  TRUE 
# 13922 16755

##---------------------------------##
## Integration  ------
##---------------------------------##

# standard normalization and variable feature selection
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000,
                                           verbose = FALSE)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)

# add proj variable
ref$proj = "Wagner2018"
seu$proj = "foxd3multiome"

# select features for downstream integration, and run PCA on each object in the list
features <- SelectIntegrationFeatures(object.list = list(ref,seu))

ref <- ScaleData(ref, features = features, verbose = FALSE)
ref <- RunPCA(ref, features = features, verbose = FALSE)

seu <- ScaleData(seu, features = features, verbose = FALSE)
seu <- RunPCA(seu, features = features, verbose = FALSE)

# find integration anchors
seu.anchors <- FindIntegrationAnchors(object.list = list(ref,seu), dims = 1:ndim)
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:ndim)


# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(seu.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = ndim, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:ndim, verbose = FALSE)

# save data
saveRDS(seu.integrated, paste0(outputdir, "integrated/multiomeall_wagner2018_cca.rds"))

# plot the original projects
p1 <- DimPlot(seu.integrated, reduction = "umap", group.by = "proj")
ggsave(filename = paste0(figdir,"seu.integrated_umap_proj.pdf"), plot = p1,width = 10, height = 7)

# plot tissue name
p1 <- DimPlot(seu.integrated, reduction = "umap", 
              group.by = "TissueName", label = TRUE, repel = TRUE)
ggsave(filename = paste0(figdir,"seu.integrated_umap_TissueName.pdf"),plot = p1, width = 15, height = 12)

# plot Wagner et al's cluster name
p1 <- DimPlot(seu.integrated, reduction = "umap", 
        group.by = "ClusterName_short", label = TRUE, repel = TRUE) 
ggsave(filename = paste0(figdir,"seu.integrated_umap_ClusterName_short.pdf"), plot = p1,width = 15, height = 12)

# plot multiome seurat clusters
p1 <- DimPlot(seu.integrated, reduction = "umap", 
        group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
ggsave(filename = paste0(figdir,"seu.integrated_umap_ClusterName_short.pdf"), plot = p1,width = 15, height = 12)

# write session info
writeLines(capture.output(sessionInfo()), "/t1-data/project/tsslab/zhu/multiome/analysis_newref/integration/sessioninfo/cca_wagner2018_multiome_sessionInfo.txt")




