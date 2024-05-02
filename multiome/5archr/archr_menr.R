# ArchR: motif enrichment analysis (Danio-code motifs)
# ZHiyuan Hu
# 4 Jan 2023
# last modified 6 Jan 2023

library(Seurat)
library(ArchR)
library(testthat)
library(parallel)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)

addArchRThreads(threads = 24) 

setwd("/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr")

ArrowFiles <- paste0("/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/","scmo_s", 1:8,".arrow")

geneAnnotation <- readRDS("Annotation/geneAnnotation_GRCz11_105.rds") 
genomeAnnotation <- readRDS("Annotation/genomeAnnotation_GRCz11_105.rds")
seqlevelsStyle(genomeAnnotation$chromSizes) <- "NCBI"
seqlevelsStyle(genomeAnnotation$blacklist) <- "NCBI"
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- "NCBI"

# load data
proj <- loadArchRProject(path = "/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/")
markersPeaks <- readRDS("/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/rds/markersPeaks.rds")

print("addMotifAnnotations")
# Danio code motifs
pwmList <- readRDS("/t1-data/project/tsslab/zhu/multiome/public/Danio_code/Regulatory_motifs_and_regulatory_site_annotations_dr11/dr11_weight_matrices_filtered_pwmList4archr.rds")
proj <- addMotifAnnotations(ArchRProj = proj, 
                            motifPWMs = pwmList,
                            name = "Motif", force = TRUE)

proj <- saveArchRProject(ArchRProj = proj)

# ### Motif enrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
saveRDS(enrichMotifs,"/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/rds/enrichMotifs.rds")

# plot heatmap
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "DanioCodeMotif-Enriched-Marker-Heatmap", width = 12, height = 6, ArchRProj = proj, addDOC = FALSE)

# write session info
writeLines(capture.output(sessionInfo()), "/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/sessionInfo/archr_menr.txt")

