#!/bin/bash
#SBATCH --job-name=footprint
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb

## Zhiyuan
## 26 Mar 2024
## last modified 27 Mar 2024

module load cuda/11.2
eval "$(conda shell.bash hook)"

conda activate chrombpnet

cluster=$1  
echo $cluster

model=${wkdir}/data/07chrombpnet_model/${cluster}/models/chrombpnet_nobias.h5
peaks=${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed
fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
chromsizes=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
fold0=${wkdir}/data/splits/fold_0.json
outdir=${wkdir}/data/18footprint_new_motifs/${cluster}/${cluster}
mkdir ${wkdir}/data/18footprint_new_motifs/${cluster}

module load bedtools/2.31.0 
bedtools intersect -v -a $peaks \
                   -b chrombpnet/data/18footprint_new_motifs/chrom_start_end_bl.bed  >  \
                   ${wkdir}/data/18footprint_new_motifs/${cluster}/${cluster}.peaks_filtered_2024march27.bed

chrombpnet footprints -m $model \
                      -r ${wkdir}/data/18footprint_new_motifs/${cluster}/${cluster}.peaks_filtered_2024march27.bed\
                      -g $fa \
                      -fl $fold0 -op $outdir \
                      -pwm_f chrombpnet/data/18footprint_new_motifs/motif_to_pwm.inhouse.tsv

conda deactivate
conda deactivate
