#!/bin/bash
#SBATCH --job-name=splitbam
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=50gb

# Zhiyuan Hu
# 17 Oct 2023
# last modified 17 Oct 2023

eval "$(conda shell.bash hook)"

conda activate split_scATAC

sample=$1
csv=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/01prepare/${sample}_cluster.csv
bam=/home/huzhiy/projects_ox/multiome/analysis_newref/cellranger_arc/scmo_${sample}/outs/atac_possorted_bam.bam

/home/huzhiy/software/pyflow-scATACseq/scripts/split_scATAC_bam_by_cluster.py \
        -prefix scmo_${sample} \
        -outdir /home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/02split_bam/scmo_${sample} $csv $bam

conda deactivate
conda deactivate