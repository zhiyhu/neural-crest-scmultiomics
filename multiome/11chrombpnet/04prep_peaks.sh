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

refdir=/home/huzhiy/projects_ox/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet

# # create chrom size file
# paste ${refdir}/chrName.txt ${refdir}/chrLength.txt > \
#       ${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes
## manually modify MT in the output to M

# oldbl=/home/huzhiy/projects_ox/ref/blacklist/danRer11_blacklist_USCSliftover_20221223.bed
# cut -c 4- $oldbl > $refdir/danRer11_blacklist_USCSliftover_20221223_noCHR.bed
# 
bedtools slop -i /home/huzhiy/projects_ox/ref/blacklist/danRer11_blacklist_USCSliftover_2023oct.bed \
              -g ${wkdir}/data/ref/GRCz11_ensembl105_foxd3.chrom.sizes \
              -b 1057 > ${wkdir}/data/ref/temp.bed
#             
# filter out blacklists    
#dNC_hoxa2b dNC_nohox_cycling NPB_nohox  mNC_arch1 mNC_arch2 mNC_head_mesenchymal mNC_hox34 mNC_nohox mNC_vagal Mutant_hox2 Mutant_nohox_12_22ss NC_trunk Pigment_gch2_high Pigment_sox6_high  dNC_hox34 dNC_nohox Mutant_hox3 Mutant_nohox_cycling Mutant_nohox_early Mutant_pigment NPB_hox2 NPB_hox3 NPB_nohox_cycling 

for cluster in merged_NPB_nohox merged_dNC_nohox mNC_head_mesenchymal mNC_arch1 Pigment_sox6_high  Mutant_nohox_early mNC_nohox Mutant_nohox_12_22ss
do
bedtools intersect -v -a ${wkdir}/data/04_1macs2/${cluster}_peaks.narrowPeak \
                   -b ${wkdir}/data/ref/temp.bed  >  \
                   ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed
done

for cluster in merged_NPB_nohox merged_dNC_nohox mNC_head_mesenchymal mNC_arch1 Pigment_sox6_high  Mutant_nohox_early mNC_nohox Mutant_nohox_12_22ss #  #
do
# sed 's/^/chr/' ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed > ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp.bed
awk -v OFS='\t' '$1 ~ /^chr[0-9]+$/ {split($1,a,"chr"); if (a[2] >= 1 && a[2] <= 25) print $0}' \
${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed > ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed

echo $cluster

cut -f1  ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed  | sort | uniq
done

# # filter chromosomes
# cd /home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/ref
# head -n 25  GRCz11_ensembl105_foxd3.chrom.sizes >  \
#             GRCz11_ensembl105_foxd3.chrom.subset.sizes
#             
            

