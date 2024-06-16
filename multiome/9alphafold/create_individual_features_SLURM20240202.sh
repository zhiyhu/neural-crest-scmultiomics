#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=array
#SBATCH --time=10:00:00

#log files:
#SBATCH -e logs/create_individual_features_%A_%a_err.txt
#SBATCH -o logs/create_individual_features_%A_%a_out.txt


#qos sets priority
#SBATCH --qos=low

#Limit the run to a single node
#SBATCH -N 1

#Adjust this depending on the node
#SBATCH --ntasks=8
#SBATCH --mem=64000

## run with sbatch --array=1-$count%100 create_individual_features_SLURM20240202.sh

eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/AlphaPulldown

fafile=/home/huzhiy/projects_ox/multiome/analysis_newref/alphafold/data/targets20240202.fasta

create_individual_features.py \
  --fasta_paths=${fafile} \
  --data_dir=/data/share_for_user/alphfold_db \
  --save_msa_files=False \
  --output_dir=/home/huzhiy/projects_ox/multiome/analysis_newref/alphafold/output/features20240128 \
  --use_precomputed_msas=False \
  --max_template_date=2050-01-01 \
  --skip_existing=True \
  --seq_index=$SLURM_ARRAY_TASK_ID