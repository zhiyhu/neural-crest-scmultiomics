#!/bin/bash
#SBATCH --job-name=predbw
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50gb

#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

## Zhiyuan
## 18 OCT 2023
## last modified 27 OCT 2023

module load cuda/11.2
eval "$(conda shell.bash hook)"

conda activate chrombpnet_test

cluster=$1  
echo $cluster

model=${wkdir}/data/07chrombpnet_model/${cluster}/models/chrombpnet_nobias.h5
peaks=${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed
fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
chromsizes=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
outdir=${wkdir}/data/08pred_bw/${cluster}

chrombpnet pred_bw -cmb $model \
                   -r $peaks  \
                   -g $fa \
                   -c $chromsizes \
                   -op $outdir
