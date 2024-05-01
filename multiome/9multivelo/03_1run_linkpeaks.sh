#!/bin/bash
#SBATCH --job-name=run_link
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=zhiyuan.hu@ndcls.ox.ac.uk     # Where to send mail  
#SBATCH -o logs/%x_%j.log 
#SBATCH -e logs/%x_%j.log 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 # specify number of CPUs to use here
#SBATCH --mem-per-cpu=5gb

## Zhiyuan
## 21 Dec 2022
## last modified 29 Dec 2023

module load R-cbrg/current

for i in {0..27000..1000}
do
sbatch -J linkp_$i /ceph/home/z/zhu/t1data/multiome/analysis_newref/multivelo2023dec/scripts/03_1signac_linkpeaks.sh $i
done