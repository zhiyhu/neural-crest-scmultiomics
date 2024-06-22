# title: "Preprocess ATAC data for multivelo"
# author: "Zhiyuan Hu"
# date: '2022-12-20'
# last modified 29 dec 2023

## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------
# knitr::opts_chunk$set(echo = TRUE)
library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomeInfoDb)

seu <- readRDS("multiome/analysis_newref/multivelo2023dec/data/02multiome/nc_multiome_forLinkPeaks.rds")
seu$stage <- as.character(seu$stage)

# HVG
hvg <- read.csv("multiome/analysis_newref/multivelo2023dec/data/wt_nohox/hvg.csv")

# normalize the gene expression data using SCTransform
DefaultAssay(seu) <- "RNA"
seu <- SCTransform(seu)

# link peaks to genes
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- 'NCBI'

# first compute the GC content for each peak
DefaultAssay(seu) <- "peaks"
seu <- RegionStats(seu, genome = BSgenome.Drerio.UCSC.danRer11)

saveRDS(seu, "multiome/analysis_newref/multivelo2023dec/data/multiome/nc_multiome_postRegionStats.rds")

## ~~~~~~~~~~~~~~~~~~~~~~~~~
## Save data --------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~

DefaultAssay(seu) <- "peaks"
seu <- LinkPeaks(
  object = seu,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = hvg$feature[hvg$highly_variable == "True"]
)

x <- Links(object = seu)
#
saveRDS(x, "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/signac_linkpeaks_hvg.rds")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Convert to 10x linkage format --------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DefaultAssay(seu) <- "peaks"

library(GenomicRanges)
library(magrittr)

# Convert the links to 10x Feature Linkage BEDPE format
# https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/analysis#tf:~:text=Feature-,Linkage,-Feature%20Linkage%20BEDPE
peak_annot <- read.delim("multiome/analysis_newref/multivelo2023dec/data/peaks/multivelo_peakanno_scplus.txt")
peak_annot$peaks <- paste0("chr",peak_annot$chrom, "-", peak_annot$start, "-", peak_annot$end)
peak_annot$peak_name <- paste0(peak_annot$gene, "_", peak_annot$peak_type)

x <- readRDS("multiome/analysis_newref/multivelo2023dec/data/wt_nohox/signac_linkpeaks_hvg.rds")

seqlevelsStyle(x) <- "UCSC"
df <- data.frame(chrom1 = 1:length(x@seqnames))
df$chrom1 <- seqnames(x) %>% as.character
df$chrom1 <- gsub(pattern = "chr", "", df$chrom1)
df$start1 <- start(x) %>% as.numeric
df$end1 <- end(x) %>% as.numeric

peaks <- x@elementMetadata$peak
df$chrom2 <- sapply(peaks, function(x) unlist(strsplit(x, "-"))[1])
df$chrom2 <- gsub(pattern = "chr", "", df$chrom2)
df$start2 <- sapply(peaks, function(x) unlist(strsplit(x, "-"))[2])
df$end2 <- sapply(peaks, function(x) unlist(strsplit(x, "-"))[3])

df$name <- x@elementMetadata$peak
df$score <- x@elementMetadata$zscore
df$strand1 <- "."
df$strand2 <- "."
df$significance <- -log10(x@elementMetadata$pvalue)
df$distance <- abs(round((as.numeric(df$start2) + as.numeric(df$end2))/2 - (as.numeric(df$start1) + as.numeric(df$end1))/2))
df$linkage_type <- "gene-peak" 

peaks_df <- paste0("chr",x@elementMetadata$peak)
# the gene-peak order in name needs to be the same as the one in the linkage type
df$name <- paste0("<",x@elementMetadata$gene,"><",peak_annot$peak_name[match(peaks_df, peak_annot$peaks)], ">")
head(df)

# write linkage bedpe file
write.table(df, "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/signac_feature_linkage_hvg.bedpe", 
            quote = FALSE, sep = "\t", row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~
## wildtype
## ~~~~~~~~~~~~~~~~~~~~~~~~~

barcodes <- read.delim( "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/cell_barcodes.txt", header = FALSE)
DefaultAssay(seu) <- "peaks"
seu[["RNA"]] <- NULL
saveRDS(seu[,barcodes$V1], "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/atac_seuobj_cistopicPeaks.rds")

library(SeuratData)
library(SeuratDisk)
SaveH5Seurat(seu[,barcodes$V1], filename = "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/atac_adata_cistopicPeaks.h5Seurat")
Convert("multiome/analysis_newref/multivelo2023dec/data/wt_nohox/atac_adata_cistopicPeaks.h5Seurat", dest = "h5ad")

## ~~~~~~~~~~~~~~~~~~~~~~~~~
## session info ------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~

# save sessionInfo for reproducibility ----
writeLines(capture.output(sessionInfo()), "multiome/analysis_newref/multivelo2023dec/sessioninfo/wtnohox_atac.txt")
