#!/bin/bash
#SBATCH --job-name=homer_annot

## Zhiyuan
## 21 Dec 2022
## last modified 28 Dec 2023

module load homer/20201202
module load bedtools/2.29.2

bed=${wkdir}/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/consensus_peak_calling/consensus_regions.bed
filteredbed=${wkdir}/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/consensus_peak_calling/consensus_regions_filtered.bed
out=${wkdir}/multiome/analysis_newref/multivelo2023dec/data/01peaks/homerannot_peaks_scplus.txt

fa=${wkdir}/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/fasta/genome.fa
gtf=${wkdir}/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf
gtfgz=${wkdir}/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/genes/genes.gtf.gz
gtfnew=${wkdir}/ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine/Danio_rerio.GRCz11.105_foxd3_mcherry_citrine.protein_coding_filtered.gtf
blist=${wkdir}/ref/blacklist/danRer11_blacklist_USCSliftover_20221223_noCHR.bed

gunzip -c  $gtfgz > $gtf
module load cellranger/7.0.1   
cellranger mkgtf $gtf $gtfnew \
                   --attribute=gene_biotype:protein_coding 

bedtools intersect -v -a $bed -b $blist  > $filteredbed

annotatePeaks.pl $filteredbed $fa -gtf $gtfnew   > $out

module purge
module load R-cbrg/current
Rscript ${wkdir}/multiome/analysis_newref/multivelo2023dec/scripts/01annot_peaks.R
