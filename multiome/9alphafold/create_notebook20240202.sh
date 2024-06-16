#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=nb
#SBATCH --time=2-00:00:00

#log files:
#SBATCH -e logs/create_notebook_%A_%a_err.txt
#SBATCH -o logs/create_notebook_%A_%a_err.txt

#Limit the run to a single node
#SBATCH -N 1

#Adjust this depending on the node
#SBATCH --ntasks=1
#SBATCH --mem=64000

eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/AlphaPulldown

cd /home/huzhiy/projects_ox/multiome/analysis_newref/alphafold/output/models20240202

module load singularity/3.8.0 
singularity exec \
    --no-home \
    --bind /home/huzhiy/projects_ox/multiome/analysis_newref/alphafold/output/models20240202:/mnt \
    /home/huzhiy/software/alpha-analysis_jax_0.4.sif \
    run_get_good_pae.sh \
    --output_dir=/mnt \
    --cutoff=100