# title: "Preprocess ATAC data for multivelo"
# author: "Zhiyuan Hu"
# date: '2022-12-20'
# last modified 28 Dec 2023

## setup
# knitr::opts_chunk$set(echo = TRUE)
library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomeInfoDb)

seu <- readRDS("multiome/analysis_newref/multivelo2023dec/data/02multiome/nc_multiome_forLinkPeaks.rds")
seu$stage <- as.character(seu$stage)

# normalize the gene expression data using SCTransform
DefaultAssay(seu) <- "RNA"
seu <- SCTransform(seu)
    
# link peaks to genes
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- 'NCBI'

# first compute the GC content for each peak
DefaultAssay(seu) <- "peaks"
seu <- RegionStats(seu, genome = BSgenome.Drerio.UCSC.danRer11)
dim(seu[["RNA"]])

saveRDS(seu, "multiome/analysis_newref/multivelo2023dec/data/02multiome/nc_multiome_postRegionStats.rds")
