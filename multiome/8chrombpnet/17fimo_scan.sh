#!/bin/bash
#SBATCH --job-name=fimo
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb

## Zhiyuan
## 4 jan 2024
## last modified 4 jan 2024

## https://github.com/kundajelab/chrombpnet/blob/a5c231fdf231bb29e9ca53d42a4c6e196f7546e8/chrombpnet/evaluation/invivo_footprints/script_new.sh#L9

module load meme_suit/5.5.5
module load bedtools/2.31.0 

sample=$1  #mNC_head_mesenchymal
echo $sample

fa=/home/huzhiy/projects_ox/ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine/Danio_rerio.GRCz11_foxd3_mcherry_citrine.dna.primary_assembly.fa
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/

oldbed=${wkdir}/data/09contribs_bw/${sample}.interpreted_regions.bed
bed=${wkdir}/data/09contribs_bw/${sample}.interpreted_regions_processed.bed
bednew=${wkdir}/data/09contribs_bw/${sample}.interpreted_regions_processed_nochr.bed
regionfa=${wkdir}/data/17fimo_in/${sample}.interpreted_regions_processed.fa

cut -f 2- $oldbed > $bed
sed 's/^chr//' $bed > $bednew

bedtools getfasta -fi $fa -bed $bednew -fo $regionfa

motif_f=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/14motif_analysis/clustered_motifs.meme
outdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/17fimo_out/${sample}
# mkdir -p $outdir
fimo --o $outdir $motif_f $regionfa

