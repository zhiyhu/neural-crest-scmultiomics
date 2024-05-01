# Create h5ad for RNA data
# Zhiyuan Hu
# 13 Jan 2023
# last modified 29 Dec 2023

library(Seurat)
library(SeuratData)
library(SeuratDisk)
seu <- LoadH5Seurat("/ceph/project/tsslab/zhu/multiome/analysis_newref/velocity/data/rds/NC_RNAvelocyto_uncorrected.h5Seurat")
dr11 <- rtracklayer::import("/ceph/project/tsslab/zhu/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf.gz")
dr11 <- dr11[dr11$gene_biotype == "protein_coding" & dr11$type == "gene",]
idx <- rownames(seu[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu[["RNA"]]@counts) %in% dr11$gene_id
table(idx)
# FALSE  TRUE 
# 2163 24773 

# remove lincRNA
seu <- seu[idx,]
seu 

# check if all colnames and rownames matches
all(colnames(seu[["spliced"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["spliced"]]@counts) == rownames(seu[["RNA"]]@counts))

all(colnames(seu[["unspliced"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["unspliced"]]@counts) == rownames(seu[["RNA"]]@counts))

all(colnames(seu[["ambiguous"]]@counts) == colnames(seu[["RNA"]]@counts))
all(rownames(seu[["ambiguous"]]@counts) == rownames(seu[["RNA"]]@counts))

# metadata containing cell type annotation
metadata <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
seu_clusters <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds")
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
# old 
# old_cell_barcodes = read.delim('/ceph/project/tsslab/zhu/multiome/analysis/multivelo2023dec/data/seurat_wnn/kept_barcodes.txt', sep='\t')
# table(colnames(seu) %in% old_cell_barcodes$x)
# table(old_cell_barcodes$x %in% colnames(seu)  )
# old_seu <-  seu[,old_cell_barcodes$x]

seu_wt <- seu[,seu$cell_type %in% c("NPB_nohox","NPB_nohox_cycling","dNC_nohox","dNC_nohox_cycling",
                                    "mNC_nohox","mNC_arch1","mNC_head_mesenchymal","Pigment_gch2_high","Pigment_sox6_high") &
                seu$genotype_new == "cit"]

dim(seu_wt)
# [1] 24773  8073
table(seu_wt$cell_type)

# wildtype cell barcode
barcode_wt <- colnames(seu_wt)
write.table(barcode_wt, "/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/cell_barcodes.txt", quote = FALSE, row.names = F, col.names = F)
# save wt rna h5ad
SaveH5Seurat(seu_wt, filename = "/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_uncorrected.h5Seurat")
Convert("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_uncorrected.h5Seurat", dest = "h5ad")
# 
# #####-------------####
# ## Corrected #####
# seu <- LoadH5Seurat("/ceph/project/tsslab/zhu/multiome/analysis_newref/velocity/data/rds/NC_RNAvelocyto.h5Seurat")
# idx <- rownames(seu[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu[["RNA"]]@counts) %in% dr11$gene_id
# table(idx)
# # FALSE  TRUE 
# # 2163 24773 
# 
# # remove lincRNA
# seu <- seu[idx,]
# seu 
# 
# # metadata containing cell type annotation
# metadata <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
# seu_clusters <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds")
# idx <- rownames(seu_clusters[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu_clusters[["RNA"]]@counts) %in% dr11$gene_id
# seu_clusters <- seu_clusters[idx,]
# 
# all(rownames(metadata) == colnames(seu))
# all(colnames(seu_clusters) == colnames(seu))
# 
# seu$cell_type <- metadata$cell_type
# seu$genotype_new <- metadata$genotype_new
# seu$stage <- as.character(seu$stage)
# # add tSNE to seu
# seu@reductions[["tsne"]]  <- seu_clusters@reductions[["tsne"]]
# 
# barcode_wt <- read.delim("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/cell_barcodes.txt", header = F) 
# SaveH5Seurat(seu[,barcode_wt$V1], filename = "/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_corrected.h5Seurat")
# Convert("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/rna_adata_corrected.h5Seurat", dest = "h5ad")

