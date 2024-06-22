#!/bin/bash
#SBATCH --job-name=splitbam

# Zhiyuan Hu
# 17 Oct 2023
# last modified 17 Oct 2023

eval "$(conda shell.bash hook)"
conda activate split_scATAC

sample=$1
csv=chrombpnet/data/01prepare/${sample}_cluster.csv
bam=cellranger_arc/scmo_${sample}/outs/atac_possorted_bam.bam

/home/huzhiy/software/pyflow-scATACseq/scripts/split_scATAC_bam_by_cluster.py \
        -prefix scmo_${sample} \
        -outdir chrombpnet/data/02split_bam/scmo_${sample} $csv $bam

conda deactivate
conda deactivate