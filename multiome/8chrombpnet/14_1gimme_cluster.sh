#!/bin/bash
#SBATCH --job-name=gimmecluster
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=2-00:00:00

## Zhiyuan
## 28 Dec 2023
## last modified 28 Dec 2023

eval "$(conda shell.bash hook)"

# conda activate /home/huzhiy/miniforge3/envs/chrombpnet
# 
# wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
# PFMall=${wkdir}/data/13motif_pfm/concatnated/all_pfm.txt
# meme2meme ${wkdir}/data/13motif_pfm/*.txt > $PFMall # merge all pfm to one file
# 
# conda deactivate

# run this ~/projects_ox/multiome/analysis_newref/chrombpnet/code/convert_meme2pfm.R
# run this ~/projects_ox/multiome/analysis_newref/chrombpnet/code/convert_meme2pfm.py

conda activate /home/huzhiy/miniforge3/envs/gimme

wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
PFMall=${wkdir}/data/13motif_pfm/concatnated/all_pfm_homerFormat_unique.txt

OUTDIR=${wkdir}/data/14motif_analysis
# the 0.999 threshold is from the chrombpnet preprint
gimme cluster -t 0.999 -N 4 $PFMall $OUTDIR

conda deactivate
conda deactivate

