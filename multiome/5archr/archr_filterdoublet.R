# Filter doublet
# ZHiyuan Hu
# 26 Dec 2022
# last modified 26 Dec 2022

library(ArchR)
setwd("/t1-data/project/tsslab/zhu/multiome/analysis_newref/archr")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Creating an ArchRProject ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("Creating an ArchRProject")

ArrowFiles <- paste0("scmo_s",1:8,".arrow")

geneAnnotation <- readRDS("Annotation/geneAnnotation_GRCz11_105.rds") 
genomeAnnotation <- readRDS("Annotation/genomeAnnotation_GRCz11_105.rds")
seqlevelsStyle(genomeAnnotation$chromSizes) <- "NCBI"
seqlevelsStyle(genomeAnnotation$blacklist) <- "NCBI"
seqlevelsStyle(BSgenome.Drerio.UCSC.danRer11) <- "NCBI"

# add annotation
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnnotation, 
  genomeAnnotation = genomeAnnotation,
  outputDirectory = "scmo_all",
  copyArrows = FALSE #This is recommened so that you maintain an unaltered copy for later usage.
)
pre_doubf <- proj@cellColData
getAvailableMatrices(proj)

# filter out the doublet
proj <- filterDoublets(ArchRProj = proj)

table(pre_doubf$Sample)
saveRDS(pre_doubf, "QualityControl/celldata_preDoubletFilter.rds")

post_doubf <- proj@cellColData
saveRDS(post_doubf, "QualityControl/celldata_postDoubletFilter.rds")
