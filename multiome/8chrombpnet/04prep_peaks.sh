#!/bin/bash
#SBATCH --job-name=prepPeaks
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=5gb

# Zhiyuan Hu
# 18 Oct 2023
# last modified 23 Oct 2023

module load bedtools/2.31.0

for cluster in merged_NPB_nohox merged_dNC_nohox mNC_head_mesenchymal mNC_arch1 Pigment_sox6_high  Mutant_nohox_early mNC_nohox Mutant_nohox_12_22ss
do
bedtools intersect -v -a ${wkdir}/data/04_1macs2/${cluster}_peaks.narrowPeak \
                   -b ${wkdir}/data/ref/temp.bed  >  \
                   ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed
done

for cluster in merged_NPB_nohox merged_dNC_nohox mNC_head_mesenchymal mNC_arch1 Pigment_sox6_high  Mutant_nohox_early mNC_nohox Mutant_nohox_12_22ss #  #
do
awk -v OFS='\t' '$1 ~ /^chr[0-9]+$/ {split($1,a,"chr"); if (a[2] >= 1 && a[2] <= 25) print $0}' \
${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed > ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed

echo $cluster

cut -f1  ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed  | sort | uniq
done