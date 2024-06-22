# ArchR
# ZHiyuan Hu
# 4 Jan2023
# last modified 26 Dec 2022

library(Seurat)
library(ArchR)
library(testthat)
library(parallel)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)

addArchRThreads(threads = 24) 

setwd("multiome/analysis_newref/archr")

ArrowFiles <- paste0("multiome/analysis_newref/archr/","scmo_s", 1:8,".arrow")

geneAnnotation <- readRDS("Annotation/geneAnnotation_GRCz11_105.rds") 
genomeAnnotation <- readRDS("Annotation/genomeAnnotation_GRCz11_105.rds")
seqlevelsStyle(genomeAnnotation$chromSizes) <- "NCBI"
seqlevelsStyle(genomeAnnotation$blacklist) <- "NCBI"
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- "NCBI"

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "scmo_all",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = F #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

### Read Seurat object
seu <- readRDS("multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx64340_clustered.rds")
cellnames <- colnames(seu)
cellnames <- gsub("_","#", cellnames)
cellnames <- gsub("s","scmo_s", cellnames)

seu$cell_type <- NA

# Neural crest
seu$cell_type[seu$seurat_clusters == 0] <- "neural crest"
seu$cell_type[seu$seurat_clusters == 2] <- "neural crest mutant"
seu$cell_type[seu$seurat_clusters == 5] <- "neural crest"
seu$cell_type[seu$seurat_clusters == 6] <- "neural crest early"
seu$cell_type[seu$seurat_clusters == 14] <- "neural crest"
seu$cell_type[seu$seurat_clusters == 15] <- "neural crest"
seu$cell_type[seu$seurat_clusters == 19] <- "neural crest"
seu$cell_type[seu$seurat_clusters == 21] <- "neural crest pigment"
seu$cell_type[seu$seurat_clusters == 24] <- "neural crest low-feature"
# posterior cluster
seu$cell_type[seu$seurat_clusters == 7] <- "notocord1"
seu$cell_type[seu$seurat_clusters == 12] <- "notocord2"
seu$cell_type[seu$seurat_clusters == 9] <- "mutant low-feature"
seu$cell_type[seu$seurat_clusters == 23] <- "hatching gland"

# tailbud
seu$cell_type[seu$seurat_clusters ==4] <- "tailbud - PSM & spinal cord"
seu$cell_type[seu$seurat_clusters ==3] <- "tailbud - PSM"
seu$cell_type[seu$seurat_clusters ==16] <- "tailbud - PSM (myotome)"
seu$cell_type[seu$seurat_clusters ==20] <- "tailbud - PSM (muscle - myl10)"

# neural 
seu$cell_type[seu$seurat_clusters == 1] <- "neural - midbrain"
seu$cell_type[seu$seurat_clusters == 8] <- "spinal cord & neurons"
seu$cell_type[seu$seurat_clusters == 13] <- "neural - diencephalon"
seu$cell_type[seu$seurat_clusters == 27] <- "differentiating neurons"

# unclassified
seu$cell_type[seu$seurat_clusters == 17] <- "unclassified1"

# other clusters
seu$cell_type[seu$seurat_clusters == 18] <- "pluripotent"

seu$cell_type[seu$seurat_clusters == 10] <- "mesoderm mixed"
seu$cell_type[seu$seurat_clusters == 11] <- "epidermal"
seu$cell_type[seu$seurat_clusters == 22] <- "endothelial"
seu$cell_type[seu$seurat_clusters == 25] <- "erythroid"

seu$cell_type[seu$seurat_clusters == 26] <- "periderm"
metadata <- seu@meta.data

# rm(seu)
table(proj$cellNames %in% cellnames)
rownames(metadata) <- cellnames

### Subset proj
proj <- proj[cellnames,]
proj <- addCellColData(ArchRProj = proj, data = metadata$cell_type,
                            cells = cellnames, name = "cell_type")

### Making Pseudo-bulk Replicates
print("addGroupCoverages")
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cell_type")

print("Saving and Loading an ArchRProject")
proj <- saveArchRProject(ArchRProj = proj)

# write session info
writeLines(capture.output(sessionInfo()), "multiome/analysis_newref/archr/sessionInfo/archr_createArchRproj.txt")



