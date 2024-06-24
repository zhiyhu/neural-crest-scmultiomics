#!/bin/sh
#SBATCH --job-name=tomtom
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log 
## Job needs 1 nodes and 1 cores per node.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 # specify number of CPUs to use here
## Job request memory
#SBATCH --mem-per-cpu=10gb

## Zhiyuan Hu
## 25 Apr 2023
## last modified 25 Apr 2023

## reference: https://github.com/bernardo-de-almeida/motif-clustering/blob/main/Motif_clustering_Drosophila.sh

#####
# Step 1: Prepare motif databases
#####

module load meme/5.4.1

# create Markov Background Model - order 3
REF=ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine
BG=GRCz11.105_foxd3_mcherry_citrine-order3.markov
fasta-get-markov -m 3 $REF $BG

# convert PWM models to MEME format
MOTIFDIR=multiome/analysis_newref/GRN_scenicplus/data/motif_with_newm/motifs_cb_format
MEMEFILE=multiome/analysis_newref/GRN_scenicplus/data/motif_with_newm/TOMTOM/all.dbs.meme
chen2meme ${MOTIFDIR}/*cb -bg $BG > $MEMEFILE

#####
# Step 2: Compute pair-wise motif similarity
#####

TOMOUT=multiome/analysis_newref/GRN_scenicplus/data/motif_with_newm/TOMTOM/tomtom.all.txt
tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	$MEMEFILE $MEMEFILE \
> $TOMOUT

# remove last lines that have details
OUTPUT=multiome/analysis_newref/GRN_scenicplus/data/motif_with_newm/TOMTOM/tomtom.all.treated.txt
head -n -4 $TOMOUT > $OUTPUT

#####
# Step 3: Hierarchically cluster motifs by similarity in R (Motif_clustering.R)
#####