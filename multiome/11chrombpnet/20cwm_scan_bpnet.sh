#!/bin/bash
#SBATCH --job-name=scan
## #SBATCH --gres=gpu:1
## #SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40gb

## Zhiyuan
## 29 OCT 2023
## last modified 29 OCT 2023

module load cuda/7.0
eval "$(conda shell.bash hook)"

cluster=mNC_head_mesenchymal

conda activate bpnet
modir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/07chrombpnet_model/${cluster}/auxiliary/interpret_subsample
contribs=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/09contribs_bw/${cluster}.profile_scores.h5

bpnet cwm-scan $modir /home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/14cwm_scan/${cluster}.csv.gz --contrib-file $contribs --num-workers 10



