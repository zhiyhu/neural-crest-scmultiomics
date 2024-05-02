# collect metadata to chrombpnet analysis
# Zhiyuan HU
# 17 Oct 2023

# refer to https://github.com/crazyhottommy/pyflow-scATACseq


# make sure the cluster_id.csv is named exactly as sample.csv

# e.g. if the sample column is sample1, it should be sample1.csv for the cluster information csv file.
# 
# cat example_cluster.csv
# name,value
# AAACGAAAGCACCATT-1,32
# AAACGAAAGTAGTCGG-1,4
# AAACGAAGTAACGGCA-1,31
# AAACTCGAGAACAGGA-1,17
# AAACTCGAGTCCAGAG-1,1
# AAACTCGCAAAGGAAG-1,41
# AAACTCGCAAGGCGTA-1,5
# AAACTCGCAGAAAGAG-1,4
# AAACTCGCAGTAGGCA-1,24
# 
# The white_list is a one column txt file with valid barcode you want to include.
# you can create it by:
#   
#   cat example_cluster.csv | cut -f1 -d, | sed '1d' > example_white_list.txt

datadir <- "/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/01prepare"
## use previous clustering results to split bam files
metadata <- readRDS("/home/huzhiy/projects_ox/multiome/analysis_newref/clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
samples <- unique(metadata$sample_id)
samples

for(i in 1:length(samples)){
  tmp <- metadata[metadata$sample_id == samples[i],]
  bc2cluster <- data.frame(name = row.names(tmp), value = tmp$cell_type)
  bc2cluster$name <- sapply(bc2cluster$name, function(x) unlist(strsplit( x,"_"))[2])

  # write cluster csv
  write.csv(bc2cluster, 
            paste0(datadir,"/",gsub("scmo_","",samples[i]),"_cluster.csv"), 
            quote = F, row.names = F) 
  # write white list csv
  write.table(bc2cluster$name, 
              paste0(datadir,"/",gsub("scmo_","",samples[i]),"_white_list.csv"), 
              quote = F, row.names = F, sep = "\n", col.names = F)
}



