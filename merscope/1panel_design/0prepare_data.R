# prepare data for Spapros to generate targeted probe set
# Zhiyuan Hu
# 16 Nov 2023
# last modified 16 Nov 2023

# prepare ss3 data
library(Seurat)
library(SeuratDisk)
library(SeuratData)
outputdir="integration/data/"

seu.integrated <- readRDS(paste0(outputdir, "integrated/NC_cca_ss3_mulitome.rds"))
ncall <- seu.integrated[,seu.integrated$method ==  "Smart-seq3"]

# Make a set of data.frames to keep track during cluster assignment.
clusters <- length(unique(ncall$seurat_clusters))
cluster.assignments <- data.frame(
  cluster=(1:clusters)-1,
  name=rep(NA, clusters),
  row.names=(1:clusters)-1
)

cluster.assignments["0","name"] <- "dNC_nohox"
cluster.assignments["5","name"] <- "NPB_nohox"
cluster.assignments["6","name"] <- "Mutant_nohox_early"
cluster.assignments["15","name"] <- "NPB_hox2"

cluster.assignments["12","name"] <- "dNC_hoxa2b"
cluster.assignments["14","name"] <- "NPB_hox3"
cluster.assignments["16","name"] <- "dNC_hox34"
cluster.assignments["17","name"] <- "NC_trunk"


cluster.assignments["8","name"] <- "Mutant_hox23"
cluster.assignments["11","name"] <- "Mutant_nohox_12_22ss"
cluster.assignments["1","name"] <- "Mutant_nohox_pigment"

cluster.assignments["10","name"] <- "mNC_vagal"
cluster.assignments["13","name"] <- "mNC_hox34"
cluster.assignments["3","name"] <- "Pigment"
cluster.assignments["2","name"] <- "mNC_nohox"
cluster.assignments["9","name"] <- "mNC_arch2"

cluster.assignments["4","name"] <- "mNC_arch1"

cluster.assignments["7","name"] <- "mNC_head_mesenchymal"

identical(ncall@assays$RNA@data, ncall@assays$RNA@counts)
ncall@assays$RNA@data <- ncall@assays$RNA@counts
DefaultAssay(ncall) <- "RNA"
dim(ncall@assays$RNA@data)

ncall@meta.data$cell_type <- NA
ncall@meta.data$cell_type <- cluster.assignments$name[match(ncall$seurat_clusters, cluster.assignments$cluster)]

SaveH5Seurat(ncall, filename ="integration/data/ss3_ncall.h5Seurat", overwrite = TRUE)
Convert("integration/data/ss3_ncall.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = TRUE)
