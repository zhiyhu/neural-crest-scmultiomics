# Create h5ad for RNA data
# Zhiyuan Hu
# 13 Jan 2023
# last modified 29 Dec 2023

library(Seurat)
library(SeuratData)
library(SeuratDisk)
seu <- LoadH5Seurat("multiome/analysis_newref/velocity/data/rds/NC_RNAvelocyto_uncorrected.h5Seurat")
dr11 <- rtracklayer::import("ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf.gz")
dr11 <- dr11[dr11$gene_biotype == "protein_coding" & dr11$type == "gene",]
idx <- rownames(seu[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu[["RNA"]]@counts) %in% dr11$gene_id

# remove lincRNA
seu <- seu[idx,]

# check if all colnames and rownames matches
all(colnames(seu[["spliced"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["spliced"]]@counts) == rownames(seu[["RNA"]]@counts))

all(colnames(seu[["unspliced"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["unspliced"]]@counts) == rownames(seu[["RNA"]]@counts))

all(colnames(seu[["ambiguous"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["ambiguous"]]@counts) == rownames(seu[["RNA"]]@counts))

# metadata containing cell type annotation
metadata <- readRDS("multiome/analysis_newref/clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
seu_clusters <- readRDS("multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds")
idx <- rownames(seu_clusters[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu_clusters[["RNA"]]@counts) %in% dr11$gene_id
seu_clusters <- seu_clusters[idx,]

all(rownames(metadata) == colnames(seu))
all(colnames(seu_clusters) == colnames(seu))

seu$cell_type <- metadata$cell_type
seu$genotype_new <- metadata$genotype_new
seu$stage <- as.character(seu$stage)
# add tSNE to seu
seu@reductions[["tsne"]]  <- seu_clusters@reductions[["tsne"]]
 
# subset for wildtype -----
seu_wt <- seu[,seu$cell_type %in% c("NPB_nohox","NPB_nohox_cycling","dNC_nohox","dNC_nohox_cycling",
                                    "mNC_nohox","mNC_arch1","mNC_head_mesenchymal","Pigment_gch2_high","Pigment_sox6_high") &
                seu$genotype_new == "cit"]

dim(seu_wt)
# [1] 24773  8073
table(seu_wt$cell_type)

# wildtype cell barcode
barcode_wt <- colnames(seu_wt)
write.table(barcode_wt, "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/cell_barcodes.txt", quote = FALSE, row.names = F, col.names = F)
# save wt rna h5ad
SaveH5Seurat(seu_wt, filename = "multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_uncorrected.h5Seurat")
Convert("multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_uncorrected.h5Seurat", dest = "h5ad")
