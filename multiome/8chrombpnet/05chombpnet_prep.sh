#!/bin/bash
#SBATCH --job-name=bpnet_prep
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem=60gb

## Zhiyuan
## 18 OCT 2023
## last modified 24 OCT 2023

# #SBATCH --gres=gpu:1
# #SBATCH --partition=gpu

eval "$(conda shell.bash hook)"

conda activate chrombpnet

wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
refdir=/home/huzhiy/projects_ox/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star

# chrombpnet prep splits
mkdir ${wkdir}/data/splits
chrombpnet prep splits -c ${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.chrom.subset.sizes \
                       -tcr chr1 chr3 chr6 \
                       -vcr chr8 chr20 -op ${wkdir}/data/splits/fold_0

## select chromosomes 1-25 in fasta
fa=${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.fa
sed '/^>/ s/>/>chr/' /home/huzhiy/projects_ox/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/fasta/genome.fa > \
    $fa 

# select chromosomes 1-25 in the peaks
for cluster in merged_NPB_nohox merged_dNC_nohox mNC_head_mesenchymal mNC_arch1 Pigment_sox6_high  Mutant_nohox_early mNC_nohox Mutant_nohox_12_22ss #  #  
do 
awk -v OFS='\t' '$1 ~ /^chr[0-9]+$/ {split($1,a,"chr"); if (a[2] >= 1 && a[2] <= 25) print $0}' \
${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist.bed > ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed

echo $cluster

# sort and remove duplicates
cut -f1  ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed  | sort | uniq
done

# chrombpnet prep nonpeaks
cluster=$1
echo $cluster
chrombpnet prep nonpeaks -g $fa \
                         -p ${wkdir}/data/04peaks/${cluster}.peaks_no_blacklist_tmp2.bed \
                         -c ${wkdir}/data/ref/GRCz11_ensembl105_foxd3.new.chrom.subset.sizes \
                         -fl ${wkdir}/data/splits/fold_0.json \
                         -br /home/huzhiy/projects_ox/ref/blacklist/danRer11_blacklist_USCSliftover_2023oct.bed \
                         -o ${wkdir}/data/05output/${cluster}


