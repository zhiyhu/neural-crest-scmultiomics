#!/usr/bin/env Rscript
# Run SoupX to filter out the ambient RNA
# Zhiyuan Hu
# 19 Aug 2023
# last updated: 29 Dec 2023

# refer to the tutorial: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html

library(SoupX)
library(DropletUtils)

##----------##
## P1       ##
##----------##

# Load data and estimate soup profile
setwd("perturb_seq/")
sc = load10X("cellranger_ccb/output/p1/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p1.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_ccb/output/p1/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()


##----------##
## P2       ##
##----------##
setwd("perturb_seq/")
# Load data and estimate soup profile
sc = load10X("cellranger/output/p2/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p2.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p2/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()

##----------##
## P3      ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_ccb/output/p3/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p3.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_ccb/output/p3/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()



##----------##
## P5       ##
##----------##

setwd("perturb_seq/")

# Load data and estimate soup profile
sc = load10X("cellranger/output/p5/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p5.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p5/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()


##----------##
## P6       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger/output/p6/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p6.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p6/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()


##----------##
## P8       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger/output/p8/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p8.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p8/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()

##----------##
## P9       ##
##----------##
setwd("perturb_seq/")

# Load data and estimate soup profile
sc = load10X("cellranger/output/p9/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p9.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p9/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()

##----------##
## P10      ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger/output/p10/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p10.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p10/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()

##----------##
## P11      ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger/output/p11/outs/")
# Estimate rho
png("perturb_seq/analysis/ambient_rna_removal/figures/soupx_rho_p11.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger/output/p11/outs/strainedCounts_soupX", out, overwrite = TRUE)
rm(sc, out)
gc()
