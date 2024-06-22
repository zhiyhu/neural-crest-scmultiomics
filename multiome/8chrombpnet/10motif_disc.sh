#!/bin/bash
#SBATCH --job-name=tfmodisco
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80gb

## Zhiyuan
## 18 OCT 2023
## last modified 10 Dec 2023

eval "$(conda shell.bash hook)"
conda activate tfmodiscolit

cluster=$1 
echo $cluster
wkdir=chrombpnet
outdir=${wkdir}/data/10tfmodisco/${cluster}
h5py=${wkdir}/data/09contribs_bw/${cluster}.counts_scores.h5

modisco motifs -i ${h5py} -n 1000000 -v -o ${outdir}/modisco_counts_results.h5
