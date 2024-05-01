#!/bin/bash
#SBATCH --job-name=footprint
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb

## Zhiyuan
## 23 OCT 2023
## last modified 23 OCT 2023


module load cuda/11.2
eval "$(conda shell.bash hook)"

conda activate chrombpnet

cluster=Mutant_nohox_12_22ss  #  mNC_nohox #  # mNC_nohox # Mutant_nohox_12_22ss
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet

model=${wkdir}/data/07chrombpnet_model/${cluster}/models/chrombpnet_nobias.h5
peaks=${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed
fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
chromsizes=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
fold0=${wkdir}/data/splits/fold_0.json
outdir=${wkdir}/data/11footprint/${cluster}/${cluster}
 
mkdir ${wkdir}/data/11footprint/${cluster}
chrombpnet footprints -m $model -r $peaks -g $fa \
                      -fl $fold0 -op $outdir \
                      -pwm_f /home/huzhiy/software/chrombpnet/chrombpnet/data/motif_to_pwm.TF.tsv

conda deactivate
conda deactivate
