# Single-cell multiomic study of neural crest

This repository contains the code for the following analysis, catelogued by the sub-projects:
* Multiome
* Smart-seq3
* Perturb-seq
* ChIP-seq
* Merscope

## 10x multiome data analysis

### 1 Preprocessing and QC

* Remove ambient RNA by SoupX: `soupX_all_samples.R`

* QC and filtering
    * `seurat_initialqc_bysample.R`
    * RNA Doublet filtering: ``
    * Intersect RNA + ATAC good-quality cells: ``
    * Merge objects: ``

### 2 Clustering and cell type annotation

* Clustering and annotating all cells: 
    * `1clustering_RNAsoupx_allcells.Rmd`
    * `2annotate_RNAsoupx_allcells.Rmd`

* Clustering and annotating NC cells: 
    * `3clustering_RNAsoupx_NC.Rmd` 
    * `4annotate_RNAsoupx_NC.Rmd`

### 3 Genotyping

* Predict genotypes in sample (75% epiboly - 4ss): `amplicon_anlaysis.R`

### 4 Integration

* Use Wagner et al (2018) data to annotate multiome all cells: `classification_wagner_multiomeall.R`

* Integrate multiome, Smart-seq3 and Wagner et al (2018) data: `cca_NC_wagner2018_multiome_SS3.R`

* Compute inter-group similarity: `NC_Wagner2018_CIDER.R`

### 5 ArchR analysis

* Create Arrow files: `archr_createArrow.R`

* Filter doublet: 

* Create ArchR project: `archr_createArchRproj.R`

* Peak calling: `archr_callpeaks.R`

* Motif enrichment analysis: `anchr_menr.R`

### 6 Velocity analysis

#### RNA velocity 

* Wild-type NC RNA velocity analysis: `scvelo_wt_uncorrected.py`

* foxd3-mutant NC RNA velocity analysis: `scvelo_mut_uncorrected.py`

#### MultiVelo 

* In `multivelo/` sub-directory

#### RegVelo

* To be added

### 7 GRN 

#### SCENIC+ analysis (`scenicplus`)



#### Regulon functional analysis (`regulon_function`)

#### GRN dynamics analysis (`GRN_dynamics`)


### 8 ChromBPNet analysis

### 9 AlphaFold analysis

### 10 Figure factory

## Smart-seq3 data analysis

### 1 Preprocessing

### 2 Clustering

### 3 Velocity


## Merscope data analysis

### 1 panel design

### 2 preprocessing and QC

### 3 clustering and mapping

## Perturb-seq data analysis

### 

## ChIP-seq data analysis


