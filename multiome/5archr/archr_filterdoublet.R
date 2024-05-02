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
# Filtering 13898 cells from ArchRProject!
#   scmo_s1 : 1705 of 13061 (13.1%)
# scmo_s2 : 3025 of 24339 (12.4%)
# scmo_s3 : 3492 of 27790 (12.6%)
# scmo_s4 : 929 of 9643 (9.6%)
# scmo_s5 : 2558 of 15996 (16%)
# scmo_s6 : 689 of 8306 (8.3%)
# scmo_s7 : 0 of 247 (0%)
# scmo_s8 : 1500 of 12248 (12.2%)

table(pre_doubf$Sample)
saveRDS(pre_doubf, "QualityControl/celldata_preDoubletFilter.rds")
# scmo_s1 scmo_s2 scmo_s3 scmo_s4 scmo_s5 scmo_s6 scmo_s7 scmo_s8 
# 13061   24339   27790    9643   15996    8306     247   12248 


post_doubf <- proj@cellColData

table(post_doubf$Sample)
saveRDS(post_doubf, "QualityControl/celldata_postDoubletFilter.rds")
# scmo_s1 scmo_s2 scmo_s3 scmo_s4 scmo_s5 scmo_s6 scmo_s7 scmo_s8 
# 11356   21314   24298    8714   13438    7617     247   10748 

