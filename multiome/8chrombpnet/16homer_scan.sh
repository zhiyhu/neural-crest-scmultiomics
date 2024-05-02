#!/bin/bash
#SBATCH --job-name=hommer_scan
# # SBATCH --gres=gpu:1
# # SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200gb

## Zhiyuan
## 4 jan 2024
## last modified 4 jan 2024

## https://github.com/kundajelab/chrombpnet/blob/a5c231fdf231bb29e9ca53d42a4c6e196f7546e8/chrombpnet/evaluation/invivo_footprints/script_new.sh#L9

eval "$(conda shell.bash hook)"

conda activate /home/huzhiy/miniforge3/envs/gimme
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/
bed=${wkdir}/data/09contribs_bw/mNC_head_mesenchymal.interpreted_regions_processed.bed
bednew=${wkdir}/data/09contribs_bw/mNC_head_mesenchymal.interpreted_regions_processed_nochr.bed

fa=/home/huzhiy/projects_ox/ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine/Danio_rerio.GRCz11_foxd3_mcherry_citrine.dna.primary_assembly.fa
outdir=${wkdir}/data/15homer_scan
mkdir -p $outdir

motif_id=${wkdir}/data/14motif_analysis/tomtom_manual_curation_chosen.txt
# motif=${wkdir}/data/14motif_analysis/clustered_chosen_motifs.txt
# run this to get homer format: ~/projects_ox/multiome/analysis_newref/chrombpnet/code/meme2homer.R
motif=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/14motif_analysis/clustered_motifs.homer.txt
sed 's/^chr//' $bed > $bednew

findMotifsGenome.pl  $bednew $fa   $outdir -find  $motif -size 200  -p 4



