# Use Wagner2018 to classify cell types of multiome cells
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

figdir="multiome/analysis_newref/integration/figures/allcells_wagner_cca/"
outputdir="multiome/analysis_newref/integration/data/"
ndim=50

##---------------------------------##
## Preprocessing both datasets ------
##---------------------------------##
seu <- readRDS("multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx64340_clustered.rds")

seu <- CreateSeuratObject(counts = seu@assays$RNA@counts,
                          meta.data = seu@meta.data)

## preprocessing Wagner 2018 data------
# https://satijalab.org/loomr/loomr_tutorial
refloom <- connect(filename = "multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018.loom", mode = "r", skip.validate = TRUE)

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
ref_metadata <- read.csv("multiome/R/ref_mapping/data/Wagner2018/public_data/WagnerScience2018_metadata.csv", row.names = 1)

# create seurat object
ref <- CreateSeuratObject(counts = full.matrix, meta.data = ref_metadata)
rm(full.matrix)
rm(ref_metadata)
refloom$close_all()

ref$ClusterName_short <- sapply(ref$ClusterName, function(x) unlist(strsplit(x, "hpf-"))[2])

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

##---------------------------------##
## Label transfer  ------
##---------------------------------##

seu.anchors <- FindTransferAnchors(reference = ref, query = seu,
                                   dims = 1:ndim, reference.reduction = "pca")

# prediction
predictions <- TransferData(anchorset = seu.anchors, refdata = ref$ClusterName_short,
                            dims = 1:ndim)
seu <- AddMetaData(seu, metadata = predictions)
seu$ClusterName_short_predicted <- seu$predicted.id

predictions <- TransferData(anchorset = seu.anchors, refdata = ref$ClusterName,
                            dims = 1:ndim)
seu <- AddMetaData(seu, metadata = predictions)
seu$ClusterName_predicted <- seu$predicted.id

predictions <- TransferData(anchorset = seu.anchors, refdata = ref$TissueName,
                            dims = 1:ndim)
seu <- AddMetaData(seu, metadata = predictions)
seu$TissueName_predicted <- seu$predicted.id

# save data
saveRDS(seu@meta.data, paste0(outputdir, "cell_type_classification/multiome64340_metadata_bywagner2018.rds"))

# sessioninfo
sessionInfo()

