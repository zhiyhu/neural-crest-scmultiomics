# ArchR: peak calling
# ZHiyuan Hu
# 4 Jan2023
# last modified 5 Jan  2023

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

proj <- loadArchRProject(path = "/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/")

### Calling peaks
print("addReproduciblePeakSet")
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "cell_type", 
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  genomeSize = 1.41e+9, #https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111374;https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1077595
  pathToMacs2 = "/t1-data/project/tsslab/zhu/.conda/envs/macs2/bin/macs2"
)
getPeakSet(proj)
proj <- addPeakMatrix(proj)

print("Saving and Loading an ArchRProject")
proj <- saveArchRProject(ArchRProj = proj)

print("getMarkerFeatures")
### Get Marker peak sets

proj@cellColData$cell_type[is.na(proj@cellColData$cell_type)] <- "N/A"
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks,"/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/rds/markersPeaks.rds")

print("addMotifAnnotations")
# CIS-BP data
motifs_list <- readRDS("/t1-data/project/tsslab/zhu/multiome/R/clustering/rds/pando/cisbp_motif.rds")
## Construction of PFM<atrixList from list of PFMatrix
library(TFBSTools)
pfmList <- do.call(PFMatrixList, motifs_list)
pwmList <- toPWM(pfmList)

names(pwmList) <- sapply(pwmList, function(x) x@name) # it needs to be names
proj <- addMotifAnnotations(ArchRProj = proj, 
                            motifPWMs = pwmList,
                            name = "Motif")
proj <- saveArchRProject(ArchRProj = proj)

### Motif enrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
saveRDS(enrichMotifs,"/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/scmo_all/rds/enrichMotifs.rds")

# write session info
writeLines(capture.output(sessionInfo()), "/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr/sessionInfo/archr_callpeaks.txt")

