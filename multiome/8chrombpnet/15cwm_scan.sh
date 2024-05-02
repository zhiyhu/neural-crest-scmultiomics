#!/bin/bash
#SBATCH --job-name=scan
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40gb
# SBATCH --dependency=afterok:37174

## Zhiyuan
## 29 OCT 2023
## last modified 29 OCT 2023

## https://github.com/kundajelab/chrombpnet/blob/a5c231fdf231bb29e9ca53d42a4c6e196f7546e8/chrombpnet/evaluation/invivo_footprints/script_new.sh#L9

module load cuda/12.2
eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/chrombpnet_test

cluster=merged_NPB_nohox # mNC_nohox
echo $cluster

# TFMDIR=/home/huzhiy/miniforge3/envs/chrombpnet_test/lib/python3.10/site-packages/chrombpnet/evaluation/invivo_footprints
TFMDIR=/home/huzhiy/software/chrombpnet_new/chrombpnet-master/chrombpnet/evaluation/invivo_footprints
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
outdir=${wkdir}/data/10tfmodisco_2024jan/${cluster}
# modir=${wkdir}/data/10tfmodisco_2024jan/${cluster}/modisco_results_allChroms_profile.hdf5
modir=${wkdir}/data/15cwm_scan/${cluster}.modisco.hdf5

# run this for indices ~/projects_ox/multiome/analysis_newref/chrombpnet/code/extract_idx_for_patterns_byIC.py
ind_file=${wkdir}/data/15cwm_scan/${cluster}_pattern_indices.txt
contribs=${wkdir}/data/09contribs_bw/${cluster}.profile_scores.h5

outdir=${wkdir}/data/15cwm_scan/${cluster}
oldbed=${wkdir}/data/09contribs_bw/${cluster}.interpreted_regions.bed
newbed=${wkdir}/data/09contribs_bw/${cluster}.interpreted_regions_processed.bed

# exit when any command fails
set -e
mkdir -p $outdir

cut -f 2- $oldbed > $newbed

indices=$(paste -sd, $ind_file) 

/home/huzhiy/miniforge3/envs/chrombpnet_test/bin/python ${TFMDIR}/tf_modiscohits.py \
   --min-ic=0 \
   --pattern-inds=${indices}  \
   --outdir=${outdir} \
	  $contribs \
	  $modir $newbed
	
	   # --pattern-inds=1,2,3,4,5 \

