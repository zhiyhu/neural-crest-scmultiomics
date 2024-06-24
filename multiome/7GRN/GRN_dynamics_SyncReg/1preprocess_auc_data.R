# preprocess auc data
# Zhiyuan Hu
# 15 jan 2024
# last modified 15 jan 2024

input_dir = 'GRN_scenicplus/ncall_2023oct_ccb/output/'

# latent time input
mv_time <- read.csv("data/wtnohox/multivelo_obs.csv")
scv_time <- read.csv("data/wtnohox/scvelo_obs.csv")
head(scv_time)

# trajectories
traj1 <- c( "mNC_head_mesenchymal", "mNC_nohox","dNC_nohox_cycling", "dNC_nohox", "NPB_nohox_cycling", "NPB_nohox")
traj2 <- c( "Pigment_gch2_high", "Pigment_sox6_high", "mNC_nohox", "dNC_nohox_cycling", "dNC_nohox", "NPB_nohox_cycling", "NPB_nohox")
traj3 <- c( "mNC_arch1", "mNC_nohox","dNC_nohox_cycling" ,  "dNC_nohox",   "NPB_nohox_cycling", "NPB_nohox" )

# keep all trajectories
idx <- mv_time$X[mv_time$cell_type %in% c(traj1, traj2, traj3)]

##~~~~~~~~~~~~~~~~~~~##
## Gene based
##~~~~~~~~~~~~~~~~~~~##
auc <- fread(paste0(input_dir, 'auc/eRegulon_AUC_gene_based.tsv.gz'))
auc[1:5,1:5]
auc_colnames <- colnames(auc)
auc$barcode <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[1])
auc$sample <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[2])
auc$cell_id <- paste0(auc$sample, "_", auc$barcode)
# auc <- data.frame(auc)
auc[1:5,1:5]
auc_cells <- auc$cell_id
auc <- auc[, 2:890]
auc_regulons <- colnames(auc)
auc <- as.data.frame(auc)
rownames(auc) <- auc_cells

df_grn <- read.csv("GRN_scenicplus/ncall_2023oct_ccb/output/cytoscape/eRegulon_metadata_filtered.csv")
df_grn 

# save the data into csv files
auc_df <- auc[match(idx, auc_cells), auc_regulons %in% df_grn$Gene_signature_name]
write.csv(x = auc_df, paste0(wkdir,"/figures/process_regulon_data/auc_gb_clean.csv"), row.names = FALSE)

auc_df <- auc[, auc_regulons %in% df_grn$Gene_signature_name]
write.csv(x = auc_df, paste0(wkdir,"/figures/process_regulon_data/eRegulon_auc_gb_filtered.csv"), row.names = FALSE)

##~~~~~~~~~~~~~~~~~~~##
## Region based
##~~~~~~~~~~~~~~~~~~~##
auc <- fread(paste0(input_dir, 'auc/eRegulon_AUC_region_based.tsv.gz'))
auc[1:5,1:5]
auc_colnames <- colnames(auc)
auc$barcode <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[1])
auc$sample <- sapply(auc$Cell, function(x) unlist(strsplit(x,"___"))[2])
auc$cell_id <- paste0(auc$sample, "_", auc$barcode)
# auc <- data.frame(auc)
auc[1:5,1:5]
auc_cells <- auc$cell_id
auc <- auc[, 2:890]
auc_regulons <- colnames(auc)
auc <- as.data.frame(auc)
rownames(auc) <- auc_cells

df_grn <- read.csv("GRN_scenicplus/ncall_2023oct_ccb/output/cytoscape/eRegulon_metadata_filtered.csv")
df_grn 

# save the data into csv files
auc_df <- auc[match(idx, auc_cells), auc_regulons %in% df_grn$Region_signature_name]
write.csv(x = auc_df, paste0(wkdir,"/figures/process_regulon_data/auc_rb_clean.csv"), row.names = FALSE)

auc_df <- auc[, auc_regulons %in% df_grn$Region_signature_name]
write.csv(x = auc_df, paste0(wkdir,"/figures/process_regulon_data/eRegulon_auc_rb_filtered.csv"), row.names = FALSE)
