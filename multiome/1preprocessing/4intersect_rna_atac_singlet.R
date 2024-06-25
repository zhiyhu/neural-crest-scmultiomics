# intersect RNA and ATAC singlet cells

dfinder <- c()
for(i in 1:4){
  sample <- c("18_22ss","12_16ss", "6_10ss","epiboly_4ss")[i]
  tmp <- readRDS(paste0("multiome/analysis_newref/preprocessing/rds/doubletfinder/metadata_soupx_", sample,".rds"))
  tmp <- tmp[,c("sample_id","DoubletFinder")]
  dfinder <- rbind(dfinder, tmp)
}

table(dfinder$DoubletFinder, dfinder$sample_id)

# scmo_s1 scmo_s2 scmo_s3 scmo_s4 scmo_s5 scmo_s6 scmo_s7 scmo_s8
# Doublet (high likelyhood)    1124    1486    2822     916    1947     430      11     835
# Doublet (low likelyhood)       80      76     203      76     115      28       1      78
# Singlet                     10671   12700   15806    8296   15094    6684     209   10950

rna_doublet <- rownames(dfinder)[dfinder$DoubletFinder != "Singlet"]

prefilter <- readRDS("multiome/analysis_newref/archr/QualityControl/celldata_preDoubletFilter.rds")
postfilter <- readRDS("multiome/analysis_newref/archr/QualityControl/celldata_postDoubletFilter.rds")

demuxlet <- prefilter
demuxlet$singlet <- FALSE
demuxlet$singlet[rownames(demuxlet) %in% rownames(postfilter)] <- TRUE

rownames(demuxlet) <- gsub("scmo_","", rownames(demuxlet))
rownames(demuxlet) <- gsub("#", "_", rownames(demuxlet))

atac_doublet <- rownames(demuxlet)[demuxlet$singlet == FALSE]

# Venn plot
x <- list(
  rna_postQC = rownames(dfinder), 
  rna_singlet= rownames(dfinder)[dfinder$DoubletFinder == "Singlet"], 
  atac_singlet = rownames(demuxlet)[demuxlet$singlet == TRUE],
  atac_postQC = rownames(demuxlet)
  
)

# devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
ggsave("QualityControl/doublet/intersection_doublet_RNA_ATAC.pdf")

shared_singlet <- rownames(dfinder)[dfinder$DoubletFinder == "Singlet"]
shared_singlet <- shared_singlet[shared_singlet %in% rownames(demuxlet)[demuxlet$singlet == TRUE]]

metadata <- dfinder[shared_singlet,]
metadata <- cbind(metadata, demuxlet[shared_singlet,])
head(metadata)
metadata <- metadata[,c(1,2,4:18)]
colnames(metadata)[17] <- "demuxlet_singlet"
saveRDS(metadata, "multiome/analysis_newref/preprocessing/rds/rna_atac_singlet/metadata.rds")

