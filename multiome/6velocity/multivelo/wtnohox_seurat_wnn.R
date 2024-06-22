# last modified 2023 dec 29
library(Seurat)
library(Signac)
library(SeuratData)
library(SeuratDisk)
setwd("multiome/analysis_newref/multivelo2023dec")
proj="wt_nohox_uncorrected"
workdir="multiome/analysis_newref/multivelo2023dec/"

# preprocess RNA-----
seu <- readRDS(paste0(workdir, "data/02multiome/nc_multiome_forLinkPeaks.rds"))

# wildtype

# keep only cells
barcodes <- read.delim(paste0(workdir, "data/",proj,"/filtered_cells.txt"), header = FALSE)
seu <- seu[,barcodes$V1]

## Only keep protein coding genes
DefaultAssay(seu) <- "RNA"
dr11 <- rtracklayer::import("ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf.gz")
dr11 <- dr11[dr11$gene_biotype == "protein_coding" & dr11$type == "gene",]
idx <- rownames(seu[["RNA"]]@counts) %in% dr11$gene_name | rownames(seu[["RNA"]]@counts) %in% dr11$gene_id
table(idx)
kept_genes <- rownames(seu[["RNA"]]@counts)[idx]
seu[["RNA"]] <- CreateAssayObject(counts = seu[["RNA"]]@counts[kept_genes,])

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, do.scale = T) # not scaled for consistency with scVelo (optionally, use SCTransform)
# seu <- SCTransform(seu)
seu <- RunPCA(seu, verbose = FALSE)
print(seu[["pca"]], dims = 1:20, nfeatures = 10)
seu <- RunUMAP(seu, dims = 1:50, reduction.name = 'umap.rna',
reduction.key = 'rnaUMAP_') # optional

DefaultAssay(seu) <- "peaks"
seu <- RunTFIDF(seu)
seu <- FindTopFeatures(seu, min.cutoff = 'q0')
seu <- RunSVD(seu)
seu <- RunUMAP(seu, reduction = 'lsi', dims = 2:50,
reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional

# find weighted nearest neighbors ----
seu <- FindMultiModalNeighbors(seu, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 30) ## changed
seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional
seu <- RunTSNE(seu, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_") # optional

# extract neighborhood graph ----
nn_idx <- seu@neighbors$weighted.nn@nn.idx
nn_dist <- seu@neighbors$weighted.nn@nn.dist
nn_cells <- seu@neighbors$weighted.nn@cell.names

# save neighborhood graph ----
setwd(paste0(workdir, "data/",proj))
write.table(nn_idx, "seurat_wnn/nn_idx.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, "seurat_wnn/nn_dist.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, "seurat_wnn/nn_cells.txt", sep = ',', row.names = F, col.names = F, quote = F)
write.table(colnames(seu), file = "seurat_wnn/kept_barcodes.txt", sep = "\t",quote = FALSE)

#Additional visualisation
# add metadata

metadata <- readRDS("multiome/analysis_newref/clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
metadata <- metadata[barcodes$V1, ]
seu$cell_type <- metadata$cell_type
seu$genotype_new <- metadata$genotype_new
p1 <- DimPlot(seu, reduction = "lsi",dims = 1:2, group.by = "cell_type")
p2 <- DimPlot(seu, reduction = "lsi",dims = 2:3, group.by = "cell_type")
p3 <- DimPlot(seu, reduction = "pca",dims = 1:2, group.by = "cell_type")
p4 <- DimPlot(seu, reduction = "pca",dims = 2:3, group.by = "cell_type")
p5 <- DimPlot(seu, reduction = "wnn.umap", group.by = "cell_type")
p6 <- DimPlot(seu, reduction = "wnn.tsne", group.by = "cell_type")
(p1 | p2 | p3)/(p4 | p5 |p6)
library(ggplot2)
ggsave(paste0(workdir, "figures/",proj,"/seurat_wnn_dimplots.pdf"),width = 18, height = 9)

