#!/usr/bin/env Rscript
# Run SoupX to filter out the ambient RNA
# Zhiyuan Hu
# 1 Sep 2022
# last updated: 20 Dec 2022

# refer to the tutorial: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html

library(SoupX)
library(DropletUtils)


##----------##
## S1       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s1/outs/")
# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s1.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("./cellranger_arc/output/scmo_s1/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()


##----------##
## S2       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s2/outs")

# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s2.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("./cellranger_arc/output/scmo_s2/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()

##----------##
## S3       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s3/outs")

# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s3.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("./cellranger_arc/output/scmo_s3/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()

##----------##
## S4       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s4/outs")
# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s4.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_arc/output/scmo_s4/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()

##----------##
## S5       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s5/outs")
# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s5.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_arc/output/scmo_s5/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()


##----------##
## S6       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s6/outs")
# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s6.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_arc/output/scmo_s6/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()


##----------##
## S8       ##
##----------##

# Load data and estimate soup profile
sc = load10X("cellranger_arc/output/scmo_s8/outs")
# Estimate rho
png("multiome/analysis_newref/ambient_rna_removal/soupx/figures/soupx_rho_s8.png")
sc = autoEstCont(sc)
dev.off()
# Clean the data
out = adjustCounts(sc)

DropletUtils:::write10xCounts("cellranger_arc/output/scmo_s8/outs/strainedCounts_soupX", out)
rm(sc, out)
gc()

