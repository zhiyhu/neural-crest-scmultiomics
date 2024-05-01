#!/bin/bash
#SBATCH --job-name=tfmodisco
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb

## Zhiyuan
## 27 OCT 2023
## last modified 26 Dec 2023

# module load cuda/11.2
eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/chrombpnet

cluster=$1
echo $1
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet

input=${wkdir}/data/10tfmodisco/${cluster}/modisco_results_profile.h5
outdir=${wkdir}/data/12modisco_report/${cluster}
motifs=/home/huzhiy/projects_ox/multiome/analysis_newref/GRN_scenicplus/data/motif_with_newm/TOMTOM/all.dbs.meme

modisco report -i $input -o $outdir  \
        -s $outdir  \
        -m $motifs \
        -n 6

conda deactivate
conda deactivate

# m A MEME file containing motifs.
