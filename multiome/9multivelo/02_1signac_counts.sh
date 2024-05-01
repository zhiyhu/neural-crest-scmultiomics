#!/bin/bash
#SBATCH --job-name=signac_counts
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=zhiyuan.hu@ndcls.ox.ac.uk     # Where to send mail  
#SBATCH -o logs/%x_%j.log 
#SBATCH -e logs/%x_%j.log 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=24:00:00

## Zhiyuan
## 21 Dec 2022
## last modified 28 Dec 2023

module load R-cbrg/current

Rscript /ceph/home/z/zhu/t1data/multiome/analysis_newref/multivelo2023dec/scripts/02_1signac_ataccounts.R
