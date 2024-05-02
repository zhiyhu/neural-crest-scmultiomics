# analysis of amplicon datta
# Predict genotypes from Miseq data
# Zhiyuan Hu
# 15 Sep 2022
# last modified 27 Dec 2022

# notes: this scripts is aimed to preprocess alevin output 

library(Seurat)
library(tximport)
library(ggplot2)

# foxd3_1_NEBNext_i701_Primer     
# foxd3_2_NEBNext_i702_Primer     
# foxd3_3_NEBNext_i703_Primer  
# foxd3mut_1_NEBNext_i704_Primer
# foxd3mut_2_NEBNext_i705_Primer
# Cherry_1_NEBNext_i706_Primer  
# Cherry_2_NEBNext_i707_Primer  

setwd("/home/z/zhu/t1data/multiome/analysis_newref/")
#------
# seurat metadata of anterior NC

seu <- readRDS("/t1-data/project/tsslab/zhu/multiome/analysis_newref/clustering/rds/seuobj_rna/seu_RNAsoupx_NC.rds")
metadata <- seu@meta.data

#---------------------------------
# read in all primers for sample 8 
mat_all <- c()

for(sample_itor in paste0("s", 1:8)){
  txi <- list()

  for(i in 1:7){
    txi[[i]] <- tximport(paste0("/home/z/zhu/t1data/multiome/runs/amplicon10x/data/alevin_output/", 
                                sample_itor, "-", i, "/alevin/"), type="alevin")$counts
    txi[[i]] <- t(as.matrix(txi[[i]]))
    colnames(txi[[i]]) <- paste0(colnames(txi[[i]]), "_primer",i)
  }
  
  barcodes <- c() # all barcodes
  for(i in 1:7){
    barcodes <- c(barcodes, rownames(txi[[i]]))
  }
  barcodes <- unique(barcodes)
  barcodes_sample <- paste0(sample_itor, "_",barcodes, "-1_",barcodes, "-1")
  
  for(i in 1:7){
    txi[[i]] <- txi[[i]][match(barcodes, rownames(txi[[i]])),]
    txi[[i]][is.na(txi[[i]])] <- 0
    rownames(txi[[i]]) <- barcodes_sample
  }
  
  # 
  mat <- c()
  for(i in 1:7){ # normalise the sum of row to 1 except for the all-0 rows
    # ))
    tmp <- txi[[i]]
    mat <- cbind(mat, tmp)
  }
  
  pdf(paste0("genotyping/figures/", sample_itor,"_heatmap.pdf"), width = 12, height = 12)
  pheatmap::pheatmap(mat, scale =  "row",show_rownames = FALSE)
  dev.off()
  
  mat_all <- rbind(mat_all, mat)
}

mat_all <- mat_all[,colSums(mat_all) >= 1]
mat_all <- mat_all[rowSums(mat_all) >= 1, ]

dim(mat_all)
# [1] 9877   12

pheatmap::pheatmap(log2(mat_all + 1), show_rownames = FALSE)

# sample dataframe
df_samples <- data.frame(barcode = rownames(mat_all), sample=NA, genotype = NA)
df_samples$sample <- substr(df_samples$barcode, 1, 2)
df_samples$genotype[df_samples$sample %in% c("s1","s3","s5","s7")] <- "cit"
df_samples$genotype[df_samples$sample %in% c("s2","s4","s6")] <- "dp"
rownames(df_samples) <- df_samples$barcode

pdf(paste0("genotyping/figures/log2mat_all_heatmap.pdf"), width = 12, height = 15)
pheatmap::pheatmap(log2(mat_all+1), show_rownames = FALSE, 
                   annotation_row = df_samples[,c("sample","genotype")])
dev.off()

# scaled mat
# mat_all_scaled <- t(apply(mat_all, 1, scale))
mat_all_scaled <- log2(mat_all+1)
rownames(mat_all_scaled) <- substr(rownames(mat_all), 1, 21)

mat_all_scaled <- mat_all_scaled[rownames(mat_all_scaled) %in% rownames(metadata),]

rownames(df_samples) <- substr(rownames(df_samples), 1, 21)
df_samples <- df_samples[rownames(df_samples) %in% rownames(metadata),]

size_n <- sum(df_samples$sample != "s8")
size_n 
# [1] 2873

## prepare the training set and the testing set
idx <- sample(1:size_n, size = round(0.6*size_n))
training_mat <- mat_all_scaled[1:size_n,][idx,]
training_res <- df_samples$genotype[1:size_n][idx]
training_res[training_res == "cit"] <- 0
training_res[training_res == "dp"] <- 1

training_set <- cbind(as.numeric(training_res), training_mat)
colnames(training_set)[1] <- "y"

testing_mat <- mat_all_scaled[1:size_n,][-idx,]
testing_res <- df_samples$genotype[1:size_n][-idx]
testing_res[testing_res == "cit"] <- 0
testing_res[testing_res == "dp"] <- 1

testing_set <- cbind(as.numeric(testing_res), testing_mat)
colnames(testing_set)[1] <- "y"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fitting the linear model -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

glm.fits <- glm(formula(paste0("y~", paste(c(0,colnames(mat_all)), collapse = "+"))),
                data = data.frame(training_set), family =binomial)
summary(glm.fits)
# 
sink("genotyping/results/glmfits_binomial.txt")
print(summary(glm.fits))
sink()  # returns output to the console

glm.probs <- predict(glm.fits, type = "response")
glm.probs[1:10]

plot(x = factor(as.character(training_set[,1])), y = glm.probs)

pdf("genotyping/figures/boxplot_linearModel_testingSet.pdf")
testing.pred <- predict(glm.fits, newdata = data.frame(testing_mat), type = "response")
plot(x = factor(as.character(testing_set[,1])), y = testing.pred)
dev.off()

## refered to https://rviews.rstudio.com/2019/03/01/some-r-packages-for-roc-curves/
library(pROC)
pdf("genotyping/figures/roc_linearModel.pdf")
pROC_obj <- roc(testing_mat[,1], testing.pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue");plot(sens.ci, type="bars")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fitting the poisson model ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glm.fits <- glm(formula(paste0("y~", paste(colnames(mat_all), collapse = "+"))),
                data = data.frame(training_set), family = poisson)
sink("genotyping/results/glmfits_poisson.txt")
summary(glm.fits)
sink()

glm.probs <- predict(glm.fits, type = "response")
glm.probs[1:10]

plot(x = factor(as.character(training_set[,1])), y = glm.probs)

pdf("genotyping/figures/boxplot_poissonModel_testingSet.pdf")
testing.pred <- predict(glm.fits, newdata = data.frame(testing_mat), type = "response")
plot(x = factor(as.character(testing_set[,1])), y = testing.pred)
dev.off()

pdf("genotyping/figures/roc_poissonModel.pdf")
pROC_obj <- roc(testing_mat[,1], testing.pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

sens.ci <- ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue");plot(sens.ci, type="bars")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# best cutoff --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# refered to https://www.projectpro.io/recipes/select-best-cutoff-point-for-problem-roc-auc-curve-r to choose the best cutoff
glm.fits <- glm(formula(paste0("y~", paste(c(0,colnames(mat_all)), collapse = "+"))),
                data = data.frame(training_set), family =binomial)
summary(glm.fits)
testing.pred <- predict(glm.fits, newdata = data.frame(testing_mat), type = "response")


library(ROCR)
ROCR_pred_test <- prediction(testing.pred, testing_res)
ROCR_perf_test <- performance(ROCR_pred_test,'tpr','fpr')
cost_perf = performance(ROCR_pred_test, "cost") 

pdf("genotyping/figures/ROCR_perf_test.pdf")
plot(ROCR_perf_test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
dev.off()

cutoff <- ROCR_pred_test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])]
cutoff
# s4_AGTCGCATCGAAGTAG-1 
# 0.5531897


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Round 1: Perform prediction in sample 8 ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s8_mat <- mat_all_scaled[(size_n+1):nrow(mat_all_scaled),]
s8.pred <- predict(glm.fits, newdata = data.frame(s8_mat), type = "response")
hist(s8.pred)

df_s8 <- data.frame(score = s8.pred)
df_s8$predicted_genotype <- NA
df_s8$predicted_genotype[df_s8$score <  cutoff ] <- "cit"
df_s8$predicted_genotype[df_s8$score >=  cutoff ] <- "dp"
rownames(df_s8) <- rownames(s8_mat)
  
dim(s8_mat)
# [1] 424  17

pdf("genotyping/figures/s8_mat_prediction.pdf")
pheatmap::pheatmap(log2(s8_mat + 1), show_rownames = FALSE, annotation_row = df_s8)
dev.off()

seu$predicted_genotype <- NA
seu$predicted_genotype[match(rownames(df_s8), colnames(seu))] <- df_s8$predicted_genotype

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Round 2: Perform further knn prediction ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seu$predicted_genotype_beforeKNN <- seu$predicted_genotype
training_set <- seu@reductions$pca@cell.embeddings[rownames(df_s8),]
testing_set <- seu@reductions$pca@cell.embeddings[!colnames(seu) %in% rownames(df_s8) & seu$sample == "4mix",]

pred <- class::knn(training_set, testing_set, df_s8$predicted_genotype)
seu$predicted_genotype[match(rownames(testing_set), colnames(seu))] <- as.character(pred)

table(seu$predicted_genotype_beforeKNN)

table(seu$predicted_genotype)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Visualisation ---------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p1 <- DimPlot(seu, group.by = "genotype", reduction = "tsne") + labs(title = "Genotypes (FACS)")
p2 <- DimPlot(seu, group.by = "predicted_genotype_beforeKNN", reduction = "tsne")+ labs(title = "Predicted genotypes (before knn)")
p3 <- DimPlot(seu, group.by = "predicted_genotype", reduction = "tsne")+ labs(title = "Predicted genotypes (after knn)")
p1 + p2 + p3
ggsave("genotyping/figures/tsne_predicted_real.pdf", width = 12, height = 4)

saveRDS(seu@meta.data[,c("predicted_genotype_beforeKNN","predicted_genotype"), drop = FALSE], "genotyping/rds/predicted_genotype.rds")

rm(mat_all, mat_all_scaled, training_mat, testing_mat)
gc()

# write session info
writeLines(capture.output(sessionInfo()), "genotyping/sessionInfo.txt")

