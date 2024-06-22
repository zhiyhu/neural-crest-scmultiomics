# title: "Preprocess ATAC data for multivelo"
# author: "Zhiyuan Hu"
# date: '2022-12-20'
# last modified 29 Dec 2023

library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomeInfoDb)
library(GenomicRanges)
library(plyranges)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Convert to 10x linkage format --------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in seq(0,27000,1000)){
  start <- i
  end <- i + 999
  tmp <- readRDS(paste0("multiome/analysis_newref/multivelo2023dec/data/03peaks/signac_linkpeaks_",start,"_",end,".rds" ))
  if(i == 0){
    x <- tmp
  } else {
    x <- bind_ranges(x,tmp)
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Convert to 10x linkage format --------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DefaultAssay(seu) <- "peaks"

library(GenomicRanges)
library(magrittr)

# Convert the links to 10x Feature Linkage BEDPE format
# https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/analysis#tf:~:text=Feature-,Linkage,-Feature%20Linkage%20BEDPE
peak_annot <- read.delim("multiome/analysis_newref/multivelo2023dec/data/01peaks/multivelo_peakanno_scplus.txt")
peak_annot$peaks <- paste0("chr",peak_annot$chrom, "-", peak_annot$start, "-", peak_annot$end)
peak_annot$peak_name <- paste0(peak_annot$gene, "_", peak_annot$peak_type)

table(paste0("chr",x@elementMetadata$peak) %in% peak_annot$peaks)
# TRUE 
# 20012 

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

# write linkage bedpe file
write.table(df, "multiome/analysis_newref/multivelo2023dec/data/03feature_linkage/signac_feature_linkage_all.bedpe", 
            quote = FALSE, sep = "\t", row.names = FALSE)
