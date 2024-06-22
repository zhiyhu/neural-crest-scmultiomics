#!/bin/bash
#SBATCH --job-name=splitbam
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=5gb

# Zhiyuan Hu
# 17 Oct 2023
# last modified 17 Oct 2023

module load samtools/1.18

# mNC_arch1 mNC_arch2 mNC_head_mesenchymal mNC_hox34 mNC_nohox mNC_vagal Mutant_hox2 Mutant_nohox_12_22ss 
# NC_trunk Pigment_gch2_high Pigment_sox6_high 
# dNC_hox34 dNC_nohox Mutant_hox3 Mutant_nohox_cycling Mutant_nohox_early Mutant_pigment
# NPB_hox2 NPB_hox3 NPB_nohox_cycling
for cluster in dNC_hoxa2b dNC_nohox_cycling NPB_nohox   #  
do
samtools merge -@ 8 \
               -r chrombpnet/data/03merged_bam/${cluster}.bam \
               chrombpnet/data/02split_bam/*/*${cluster}.bam
done