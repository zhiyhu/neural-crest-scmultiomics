#!/bin/sh
#SBATCH --job-name=tomtom  
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log 
## Job needs 1 nodes and 1 cores per node.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
## Job request memory
#SBATCH --mem-per-cpu=10gb

## Zhiyuan Hu
## 30 Dec 2023
## last modified 30 Dec 2023

# reference chrombpnet biorxiv

eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/gimme

wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref
MEMEFILE=${wkdir}/chrombpnet/data/14motif_analysis/clustered_motifs.meme
DBMOTIF=${wkdir}/GRN_scenicplus/data/motif_with_newm/TOMTOM/all.dbs.meme
TOMOUT=${wkdir}/chrombpnet/data/14motif_analysis/tomtom/clustered_motifs_vs_daniocode_with_newm.txt

cd ${wkdir}/chrombpnet/data/14motif_analysis/tomtom/
tomtom \
	-no-ssc -oc . -verbosity 4 -text -min-overlap 5 -dist pearson -evalue -thresh 10.0 \
	$MEMEFILE $DBMOTIF \
  > $TOMOUT

