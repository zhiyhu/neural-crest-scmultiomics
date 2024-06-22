# title: "Preprocess ATAC counts matrix for multivelo"
# author: "Zhiyuan Hu"
# date: '2022-12-20'
# last modified 28 Dec 2023

## ----setup, include=FALSE
# knitr::opts_chunk$set(echo = TRUE)
library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomeInfoDb)

# Load and preprocess data

print("Load and preprocess data")
seu <- readRDS("multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds") # Clustered and annotated dataset for anterior NC (MUT+WT)
seu_atac <- readRDS("multiome/analysis_newref/preprocessing/rds/rna_atac_singlet/seuobj_atac64340.rds")


seu_atac <- seu_atac[,colnames(seu)]
DefaultAssay(seu) <- "RNA"
seu[["ATAC"]] <- seu_atac[["ATAC"]]

## Add gene annotation

# Import the gtf gene annotations
library(ensembldb)
library(AnnotationHub)

dr11 <- rtracklayer::import("ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf.gz")
genome(dr11) <- "GRCz11"

# convert to UCSC style
seqnames(dr11)
dr11 <- dr11[seqnames(dr11) %in% c(1:25,"MT")] ##to avoid the error: cannot switch some of GRCz11's seqlevels from NCBI to UCSC style
seqlevelsStyle(dr11) <- "NCBI" # change seqnames from e.g. chr1 to 1
seqlevelsStyle(seu[["ATAC"]]@ranges) <- "NCBI"
# set gene annotations
DefaultAssay(seu) <- "ATAC"
Annotation(seu) <- dr11
rm(dr11); gc()

# Load cisTopic-called peaks

# refer to: https://stuartlab.org/signac/articles/merging.html
# read in peak sets
print("Load cisTopic-called peaks")
scplus_peaks <- read.table(
  file = "multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/consensus_peak_calling/consensus_regions.bed",
  col.names = c("chr", "start", "end","cell_type","1","2"), row.names = NULL
)
colnames(scplus_peaks) <- c("chr","start", "end","cell_type","X1","X2")
scplus_peaks <- makeGRangesFromDataFrame(scplus_peaks)

# Quantify counts in each peak (takes ~5h)

print("Quantify counts in each peak")

# to fix the change in directory name# on 28 dec 2023
tmp <- Fragments(seu)
for(i in 1:length(tmp)){
  tmp[[i]]@path <- gsub("t1-data","ceph",tmp[[i]]@path)
}
seu <- SetAssayData(seu, slot = "fragments", new.data = tmp)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(seu),
  features = scplus_peaks,
  cells = colnames(seu)
)
saveRDS(macs2_counts, "multiome/analysis_newref/multivelo2023dec/data/02atac/counts_cistopicPeaks.rds")

chrom_assay <- CreateChromatinAssay(
  counts = macs2_counts,
  ranges = scplus_peaks,  
  fragments = seu[["ATAC"]]@fragments,
  annotation = seu[["ATAC"]]@annotation
)

DefaultAssay(seu) <- "RNA"
seu[["peaks"]] <- chrom_assay
seu[["ATAC"]]  <- NULL
saveRDS(seu, "multiome/analysis_newref/multivelo2023dec/data/02multiome/nc_multiome_forLinkPeaks.rds")

sessionInfo()
