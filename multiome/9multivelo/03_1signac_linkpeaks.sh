#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=zhiyuan.hu@ndcls.ox.ac.uk     # Where to send mail  
#SBATCH -o logs/%x_%j.log 
#SBATCH -e logs/%x_%j.log 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=10gb

## Zhiyuan
## 21 Dec 2022
## last modified 29 dec 2023

module load R-cbrg/current

echo $1
Rscript /ceph/home/z/zhu/t1data/multiome/analysis_newref/multivelo2023dec/scripts/03_1signac_linkpeaks.R $1

