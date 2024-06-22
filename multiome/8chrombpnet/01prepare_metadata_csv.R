# collect metadata to chrombpnet analysis
# Zhiyuan HU
# 17 Oct 2023

# refer to https://github.com/crazyhottommy/pyflow-scATACseq

datadir <- "chrombpnet/data/01prepare"
## use previous clustering results to split bam files
metadata <- readRDS("clustering/rds/metadata/seu_RNAsoupx_NC_metadata.rds")
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



