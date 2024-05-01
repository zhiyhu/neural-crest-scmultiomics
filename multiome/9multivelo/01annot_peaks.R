# Prepare peaks annotation for MultiVelo
# Zhiyuan Hu
# 21 Dec 2022
# last modified 28 Dec 2023

# https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/annotation
library(GenomicRanges)

# Read hommer annotation
# /ceph/project/tsslab/zhu/multiome/R/multivelo2023dec/scripts/homer_annotatepeaks.sh
peak_annot <- read.delim("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/01peaks/homerannot_peaks_scplus.txt")

print("Load cisTopic-called peaks")
scplus_peaks <- read.table(
  file = "/ceph/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/consensus_peak_calling/consensus_regions_filtered.bed",
  col.names = c("chr", "start", "end","cell_type","1","2"), row.names = NULL
)
colnames(scplus_peaks) <- c("chr","start", "end","cell_type","X1","X2")
scplus_peaks <- makeGRangesFromDataFrame(scplus_peaks)

# 1	chrom	Contig that contains the peak
# 2	start	Peak start location
# 3	end	Peak end location
# 4	gene	Gene symbol based on the gene annotation in the reference.
# 5	distance	Distance of peak from TSS of gene. Positive distance means the start of the peak is downstream of the position of the TSS, whereas negative distance means the end of the peak is upstream of the TSS. Zero distance means the peak overlaps with the TSS or the peak overlaps with the transcript body of the gene.
# 6	peak_type	Can be "promoter", "distal" or "intergenic".

df <- data.frame(chrom = peak_annot$Chr,
                 start = peak_annot$Start -1,
                 end = peak_annot$End,
                 gene = peak_annot$Gene.Name,
                 distance = peak_annot$Distance.to.TSS,
                 peak_type = peak_annot$Annotation)

# If a peak overlaps the body of a transcript, and it is not a promoter nor a distal peak of the gene, it will be annotated as a distal peak of that gene with distance set as zero.

# The annotation procedure is as follows:
#   
# If a peak overlaps with promoter region (-1000 bp, +100 bp) of any transcription start site (TSS), it is annotated as a promoter peak of the gene.
# If a peak is within 200 kb of the closest TSS, and if it is not a promoter peak of the gene of the closest TSS, it will be annotated as a distal peak of that gene.
# If a peak overlaps the body of a transcript, and it is not a promoter nor a distal peak of the gene, it will be annotated as a distal peak of that gene with distance set as zero.
# If a peak has not been mapped to any gene at the step, it will be annotated as an intergenic peak without a gene symbol assigned.

df$distance[grepl("exon", df$peak_type)] <- 0
df$peak_type[grepl("exon", df$peak_type)] <- "distal"
df$distance[grepl("intron", df$peak_type)] <- 0
df$peak_type[grepl("intron", df$peak_type)] <- "distal"
df$distance[grepl("TTS", df$peak_type)] <- 0
df$peak_type[grepl("TTS", df$peak_type)] <- "distal"
# If a peak overlaps with promoter region (-1000 bp, +100 bp) of any transcription start site (TSS), it is annotated as a promoter peak of the gene.
df$peak_type[grepl("promoter", df$peak_type)] <- "promoter"

df$peak_type[grepl("Intergenic", df$peak_type)  & abs( df$distance) <= 200000] <- "distal"
# If a peak has not been mapped to any gene at the step, it will be annotated as an intergenic peak without a gene symbol assigned.
df$peak_type[grepl("Intergenic", df$peak_type) ] <- "intergenic"

table(df$peak_type)
# distal intergenic   promoter 
# 309667       4204      23974 

range(df$distance[df$peak_type == "promoter"])
# [1] -1000   981
df <- df[!is.na(df$peak_type),]

write.table(df, sep = "\t", quote = FALSE, file = "/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/01peaks/multivelo_peakanno_scplus.txt", row.names = FALSE)


scplus_peaks_ranges <- paste(as.character(scplus_peaks@seqnames), scplus_peaks@ranges@start-1, scplus_peaks@ranges@start + 499, sep="-")
peak_annot_ranges <- paste(df$chrom, df$start, df$end, sep = "-")
table(peak_annot_ranges %in% scplus_peaks_ranges)
# FALSE   TRUE 
# 335414   2431

# library(GenomicRanges)
# grtest.unstrand <- GRanges(
#   seqnames = "6",
#   ranges=IRanges(
#     start=c(32090000), # note the difference
#     width=110000))
# 
# GenomicRanges::intersect(scplus_peaks, grtest.unstrand, )[1:5,]

df[df$chrom == 6 & df$start >= 32090000 & df$start < 32200000,]
