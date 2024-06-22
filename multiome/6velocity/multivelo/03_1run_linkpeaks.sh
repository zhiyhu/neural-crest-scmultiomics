#!/bin/bash
#SBATCH --job-name=run_link

## Zhiyuan
## 21 Dec 2022
## last modified 29 Dec 2023

module load R-cbrg/current

for i in {0..27000..1000}
do
sbatch -J linkp_$i multiome/analysis_newref/multivelo2023dec/scripts/03_1signac_linkpeaks.sh $i
done