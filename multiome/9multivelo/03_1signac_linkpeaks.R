# title: "Preprocess ATAC data for multivelo"
# author: "Zhiyuan Hu"
# date: '2022-12-20'
# last modified 29 Dec 2023

## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------
# knitr::opts_chunk$set(echo = TRUE)
library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomeInfoDb)

seu_atac <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/02multiome/nc_multiome_postRegionStats.rds")
seu_atac

args = commandArgs(trailingOnly=TRUE)
start = as.numeric(args[1])
end = start + 999
end = min(end, nrow(seu_atac))
# link peaks to genes
DefaultAssay(seu_atac) <- "RNA"
seu_atac <- LinkPeaks(
  object = seu_atac,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = rownames(seu_atac)[start:end]
)

DefaultAssay(seu_atac) <- "peaks"
x <- Links(object = seu_atac)

print("linking done.. save links")
saveRDS(x, paste0("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/03peaks/signac_linkpeaks_",start,"_",end,".rds" ))

sessionInfo()