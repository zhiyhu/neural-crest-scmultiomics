# Integrate Multiome data, NC only, and SS3 data with Wagner 2018
# Zhiyuan Hu
# 9 Sep 2023
# last modified 11 Sep 2023
## refer to https://satijalab.org/seurat/articles/integration_large_datasets.html

## prepare multiome data

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(loomR)
library(SeuratWrappers)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(patchwork)

stage_cols <- RColorBrewer::brewer.pal("Spectral", n = 11)
genotype_cols <- RColorBrewer::brewer.pal("Set2", n = 8)[c(1,2,3,8)]
sample_cols <- RColorBrewer::brewer.pal("Accent", n = 8)
doublet_cols <- RColorBrewer::brewer.pal("Paired", n = 8)[c(6,8,2)]
col.ls <- unique(c(ArchR::ArchRPalettes$stallion, ArchR::ArchRPalettes$bear)) # archR color palettes
names(col.ls) <- NULL

figdir="/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/figures/NC_3datasets_cca/"
outputdir="/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/data/"
ndim <- 50

ncmo <- readRDS("../../clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds") # multiome
ncss <- readRDS("../data/ss3/nc_anterior_seuobj_batches1to8.rds") # smartseq3

ncmo$method <- "10xmultiome"
ncss$method <- "Smart-seq3"

dim(ncmo)
# [1] 27599 16550

dim(ncss)
# [1] 27599  1181

table(rownames(ncmo) %in% rownames(ncss))
#  TRUE 
# 27599

rownames(ncmo)[!rownames(ncmo) %in% rownames(ncss)]
# [1] "citrine" "mCherry"

# make meta data columns consistent
colnames(ncss@meta.data) <- tolower(colnames(ncss@meta.data))
colnames(ncmo@meta.data) <- tolower(colnames(ncmo@meta.data))
ncss$genotype <- tolower(ncss$genotype)

# merge two objects
nc_all <- merge(ncmo, y = ncss)
# nc_all@meta.data <- nc_all@meta.data[,1:20]
ref$method <- "Wagner2018"
nc_all <- merge(nc_all, y = ref)

## Read Wagner data

## preprocessing Wagner 2018 data------
# https://satijalab.org/loomr/loomr_tutorial
refloom <- connect(filename = "/ceph/project/tsslab/zhu/multiome/analysis/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.loom", mode = "r", skip.validate = TRUE)
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
ref_metadata <- read.csv("/ceph/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018_metadata.csv", row.names = 1)

# create seurat object
ref <- CreateSeuratObject(counts = full.matrix, meta.data = ref_metadata)
rm(full.matrix)
rm(ref_metadata)
refloom$close_all()

ref$ClusterName_short <- sapply(ref$ClusterName, function(x) unlist(strsplit(x, "hpf-"))[2])
# SaveH5Seurat(ref, filename = "/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.processed.h5Seurat")
# Convert("/t1-data/project/tsslab/zhu/multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.processed.h5Seurat", dest = "h5ad")

## gene overlap
table(rownames(ref) %in% rownames(nc_all))
# FALSE  TRUE 
# 13970 16707 

## Integration

# split the dataset into a list of two seurat objects (stim and CTRL)
seu_lst <- SplitObject(nc_all, split.by = "method")

# normalize and identify variable features for each dataset independently
seu_lst <- lapply(X = seu_lst, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu_lst)

# find integration anchors
seu.anchors <- FindIntegrationAnchors(object.list = seu_lst, dims = 1:ndim)
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:ndim)


DefaultAssay(seu.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 50, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:30)
seu.integrated <- FindNeighbors(seu.integrated, reduction = "pca", dims = 1:30)
seu.integrated <- FindClusters(seu.integrated, resolution = 0.5)

# save data
saveRDS(seu.integrated, paste0(outputdir, "integrated/NC_cca_3datasets.rds"))

# plot by methods
DimPlot(seu.integrated, reduction = "umap", group.by = "method",pt.size = 0.5, cols = alpha(genotype_cols, 0.3)) # split.by = "method"
ggsave("../figures/NC_3datasets_cca/umap_by_methods.pdf", width = 9, height = 6)

# plot by methods and splitted by methods
DimPlot(seu.integrated, reduction = "umap", group.by = "method", pt.size = 0.5, cols = alpha(genotype_cols, 0.3), split.by = "method") #
ggsave("../figures/NC_3datasets_cca/umap_by_methods_splitted_by_methods.pdf", width = 18, height = 5)

# # plot the original projects
# p1 <- DimPlot(seu.integrated, reduction = "umap", group.by = "proj")
# ggsave(filename = paste0(figdir,"seu.integrated_umap_proj.pdf"), plot = p1,width = 10, height = 7)
# 
# # plot tissue name
# p1 <- DimPlot(seu.integrated, reduction = "umap", 
#               group.by = "TissueName", label = TRUE, repel = TRUE)
# ggsave(filename = paste0(figdir,"seu.integrated_umap_TissueName.pdf"),plot = p1, width = 15, height = 12)
# 
# # plot Wagner et al's cluster name
# p1 <- DimPlot(seu.integrated, reduction = "umap", 
#         group.by = "ClusterName_short", label = TRUE, repel = TRUE) 
# ggsave(filename = paste0(figdir,"seu.integrated_umap_ClusterName_short.pdf"), plot = p1,width = 15, height = 12)
# 
# # plot multiome seurat clusters
# p1 <- DimPlot(seu.integrated, reduction = "umap", 
#         group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
# ggsave(filename = paste0(figdir,"seu.integrated_umap_ClusterName_short.pdf"), plot = p1,width = 15, height = 12)
# 
# # write session info
# writeLines(capture.output(sessionInfo()), "/t1-data/project/tsslab/zhu/multiome/analysis_newref/integration/sessioninfo/cca_wagner2018_multiome_sessionInfo.txt")




