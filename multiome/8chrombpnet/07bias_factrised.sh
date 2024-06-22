#!/bin/bash
#SBATCH --job-name=biasfac
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb

## Zhiyuan
## 18 OCT 2023
## last modified 25 OCT 2023

module load cuda/11.2
eval "$(conda shell.bash hook)"

refdir=ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star

conda activate chrombpnet_test

cluster=$1 

fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
chromsizes=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
peaks=${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed
neg=${wkdir}/data/05output/${cluster}_negatives.bed
fold0=${wkdir}/data/splits/fold_0.json
outdir=${wkdir}/data/07chrombpnet_model/${cluster}/
bmod=chrombpnet/data/06bias_model/${cluster}_gputest/models/${cluster}_bias.h5
new_bam=chrombpnet/data/03merged_bam/${cluster}_new.bam

chrombpnet pipeline \
        -ibam $new_bam \
        -d "ATAC" \
        -g $fa \
        -c $chromsizes \
        -p $peaks \
        -n $neg \
        -fl $fold0 \
        -b $bmod \
        -o $outdir
