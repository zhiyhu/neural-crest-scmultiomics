# Initial QC of multiome data (filtered by soupx)
# Zhiyuan HU
# 22 Dec 2022
# last modified 22 Dec 2022

library(dplyr)
library(Seurat)
library(Signac)
library(patchwork)
library(scales)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE) # should be s1,s2 etc.
id <- paste0("scmo_", args[1])

seu <- readRDS(paste0("~/t1data/multiome/analysis_newref/preprocessing/rds/scmo_seu_obj/scmo_",args[1],"_multiome_soupxFiltered.rds"))
sample_info <- read.delim("/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/code/1seu_initial_qc/sample_table.txt", sep = "\t")
seu$sample_id <- id
seu$sample <- sample_info$sample_name[sample_info$sample_id == id]
seu$stage <- sample_info$stage[sample_info$sample_id == id]
seu$genotype <- sample_info$genotype[sample_info$sample_id == id]

# rename cells
RenameCells(seu, add.cell.id =  paste0(args[1], "_", colnames(seu)))

cat(paste0("Number of cells pre-filtering in sample ",print(id), ": ", ncol(seu)),
    file="/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/reports/seurat_initialqc_statistics.txt",
    sep="\n",append=TRUE)

p1 <- VlnPlot(
  object = seu,
  features = c("nCount_RNA", "nCount_ATAC", "nFeature_RNA", 
               "nFeature_ATAC", "nucleosome_signal", "percent.mt"),
  ncol = 6,
  pt.size = 0
)
ggsave(plot = p1, filename = paste0("/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/figures/seurat_depthQC/depth_prefilter_", id, ".pdf"), width = 18, height = 7)

# filter out low quality cells
seu <- subset(
  x = seu,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 15000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 & 
    nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    nucleosome_signal < 2 &
    percent.mt < 20
  
)

cat(paste0("Number of cells post-filtering in sample ",print(id), ": ", ncol(seu)),
    file="/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/reports/seurat_initialqc_statistics.txt",
    sep="\n",append=TRUE)

p1 <- VlnPlot(
  object = seu,
  features = c("nCount_RNA", "nCount_ATAC", "nFeature_RNA", 
               "nFeature_ATAC", "nucleosome_signal", "percent.mt"),
  ncol = 6,
  pt.size = 0
)
ggsave(plot = p1, filename = paste0("/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/figures/seurat_depthQC/depth_postfilter_", id, ".pdf"), width = 18, height = 7)

plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "sample")

ggsave(plot = (plot1 + plot2), filename = paste0("/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/figures/seurat_depthQC/featurescatter_postfilter_", id, ".pdf"), width = 10, height = 5)

# save data
saveRDS(seu, paste0("/t1-data/project/tsslab/zhu/multiome/analysis_newref/preprocessing/rds/intermediate_seuobj/post_initialqc_",id,".rds"))

sessionInfo()
