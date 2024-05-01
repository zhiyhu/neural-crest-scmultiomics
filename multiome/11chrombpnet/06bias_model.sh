#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb

## Zhiyuan
## 18 OCT 2023
## last modified 23 OCT 2023

module load cuda/11.2
eval "$(conda shell.bash hook)"

wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
refdir=/home/huzhiy/projects_ox/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star

conda activate chrombpnet_test

cluster=$1  #  mNC_nohox # Mutant_nohox_12_22ss # mNC_nohox # 
echo $cluster

wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet

bam=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/03merged_bam/${cluster}.bam
new_bam=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/03merged_bam/${cluster}_new.bam
fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
chromsizes=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
peaks=${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed
neg=${wkdir}/data/05output/${cluster}_negatives.bed
fold0=${wkdir}/data/splits/fold_0.json
outdir=${wkdir}/data/06bias_model/${cluster}_gputest
sorted_bam=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/03merged_bam/${cluster}_sorted.bam

# module load samtools
# samtools sort $new_bam -o $sorted_bam


chrombpnet bias pipeline \
        -ibam $new_bam \
        -d "ATAC" \
        -g $fa \
        -c $chromsizes \
        -p $peaks \
        -n $neg \
        -fl $fold0 \
        -b 0.5 \
        -o $outdir \
        -fp ${cluster}
        

