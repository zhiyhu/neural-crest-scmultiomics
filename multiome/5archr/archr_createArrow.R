#!/usr/bin/Rscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run archR qc
# Zhiyuan Hu
# 23 Dec 2022
# last modified 25 Dec 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# refer to this protocol: https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up -----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ArchR)
library(testthat)
library(parallel)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)

addArchRThreads(threads = 24) 
samples=1:8

setwd('multiome/analysis_newref/archr/')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating Arrow Files--------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("Creating Arrow Files")

# this is pre-run
inputFiles <- paste0("multiome/analysis_newref/cellranger_arc/output/scmo_s", samples, "/outs/atac_fragments.tsv.gz")
names(inputFiles) <- paste0("scmo_s", samples)
# see https://github.com/GreenleafLab/ArchR/issues/302
# input files need to be in gz format
# reformatFragmentFiles(inputFiles)

geneAnnotation <- readRDS("Annotation/geneAnnotation_GRCz11_105.rds") 
genomeAnnotation <- readRDS("Annotation/genomeAnnotation_GRCz11_105.rds")
seqlevelsStyle(genomeAnnotation$chromSizes) <- "NCBI"
seqlevelsStyle(genomeAnnotation$blacklist) <- "NCBI"
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- "NCBI"

# this is pre-run
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 2500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# Inferring Doublets-------

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE
)
saveRDS(doubScores,"QualityControl/doubScores.rds")

# sessioninfo-------
sessionInfo()
