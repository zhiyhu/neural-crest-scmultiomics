#!/bin/bash

## Zhiyuan
## 21 Dec 2022
## last modified 29 dec 2023

module load R-cbrg/current

echo $1
Rscript multiome/analysis_newref/multivelo2023dec/scripts/03_1signac_linkpeaks.R $1

