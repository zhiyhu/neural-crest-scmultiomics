#!/bin/bash
#SBATCH --job-name=tfmodisco
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb

## Zhiyuan
## 27 dec 2023
## last modified 27 dec 2023

eval "$(conda shell.bash hook)"
conda activate tfmodiscolit

cluster=$1
echo $1

input=${wkdir}/data/10tfmodisco/${cluster}/modisco_results_profile.h5
outdir=${wkdir}/data/13motif_pfm/${cluster}.txt

modisco meme -i $input \
        -o $outdir  \
        -t PFM

conda deactivate
conda deactivate

