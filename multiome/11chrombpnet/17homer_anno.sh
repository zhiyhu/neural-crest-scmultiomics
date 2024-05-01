#!/bin/bash
#SBATCH --job-name=homer_annot
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=50gb
#SBATCH --time=04:00:00

## Zhiyuan
## 21 Dec 2022
## last modified 10 jan 2024

eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/gimme

wkdir=/home/huzhiy/projects_ox
bed=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/17fimo_out/mNC_nohox/fimo_out_bed.tsv

out=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet/data/17fimo_out/mNC_nohox/fimo_out.homerannot_peaks.txt

fa=${wkdir}/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/fasta/genome.fa

gtfnew=${wkdir}/ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine/Danio_rerio.GRCz11.105_foxd3_mcherry_citrine.protein_coding_filtered.gtf

annotatePeaks.pl $bed $fa -gtf $gtfnew   > $out

conda deactivate
conda deactivate