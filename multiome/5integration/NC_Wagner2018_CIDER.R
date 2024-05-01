## Apply CIDER between multiome data and Wagner 2018 data
## Zhiyuan Hu
## 28 March 2023
## last modified 1 Apr 2023

library(CIDER)
library(Seurat)
library(parallel)
library(cowplot)

# read in integrated data
seu.integrated <- readRDS("/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/data/integrated/multiomeNC_wagner2018_cca.rds")
seu.integrated$initial_cluster <- NA
seu.integrated$initial_cluster <- seu.integrated$ClusterName_short
seu.integrated$initial_cluster[seu.integrated$proj =="foxd3multiome"] <- seu.integrated$cell_type[seu.integrated$proj =="foxd3multiome"]
# remove not legid symbols
seu.integrated$initial_cluster <- gsub(" ", "_", seu.integrated$initial_cluster)
seu.integrated$initial_cluster <- gsub("-", "_", seu.integrated$initial_cluster)
seu.integrated$initial_cluster <- gsub("___", "_", seu.integrated$initial_cluster)
seu.integrated$initial_cluster <- gsub("__", "_", seu.integrated$initial_cluster)
seu.integrated$initial_cluster[grep("differentiating_neurons",seu.integrated$initial_cluster)] <- "differentiating_neurons"
seu.integrated$initial_cluster[grep("epidermal",seu.integrated$initial_cluster)] <- "epidermal"
seu.integrated$initial_cluster[grep("pharyngeal_arch",seu.integrated$initial_cluster)] <- "pharyngeal_arch"
seu.integrated$initial_cluster[grep("diencephalon",seu.integrated$initial_cluster)] <- "diencephalon"
seu.integrated$initial_cluster[grep("pronephric_duct",seu.integrated$initial_cluster)] <- "pronephric_duct"
seu.integrated$initial_cluster[grep("mesoderm",seu.integrated$initial_cluster)] <- "mesoderm"
seu.integrated$initial_cluster[grep("heart",seu.integrated$initial_cluster)] <- "heart"
seu.integrated$initial_cluster[grep("muscle",seu.integrated$initial_cluster)] <- "muscle"
seu.integrated$initial_cluster[grep("lateral_line",seu.integrated$initial_cluster)] <- "lateral_line"
seu.integrated <- seu.integrated[,!grepl("apoptotic",seu.integrated$initial_cluster)]

table(seu.integrated$initial_cluster)

# run IDER
ider <- getIDEr(seu.integrated, 
                group.by.var = "initial_cluster",
                batch.by.var = "proj",
                downsampling.size = 35, downsampling.include = FALSE,
                use.parallel = TRUE, verbose = TRUE, n.cores = 6)
saveRDS(object = ider, 
        file = "/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/data/integrated/ider.rds")
# Error in seq_len(i - 1) : 
#   argument must be coercible to non-negative integer

# plot hirerchical clustering
hc <- hclust(as.dist(1-(ider[[1]] + t(ider[[1]])))/2)

pdf("/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/figures/NC_wagner_cca/cider_hclust.pdf",width = 12, height = 6)
plot(hc)
# plot(as.dendrogram(hc), horiz = F)
dev.off()

## plot the heatmap
library(viridis)
idx1 <- unique(ider[[3]]$g1[ider[[3]]$b1 == "Wagner2018"])
idx2 <- unique(ider[[3]]$g2[ider[[3]]$b2 == "foxd3multiome"])
pheatmap::pheatmap(
  ider[[1]][idx1, idx2],
  color = inferno(10),
  border_color = NA,
  display_numbers = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  width = 14,
  height = 24,
  filename =  "/ceph/project/tsslab/zhu/multiome/analysis_newref/integration/figures/NC_wagner_cca/cider_heatmap.pdf"
)
