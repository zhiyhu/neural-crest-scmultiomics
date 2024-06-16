# Single-cell multiomic study of neural crest

This repository contains the code for the following analysis, catelogued by the sub-projects:
* Multiome
* Smart-seq3
* Perturb-seq
* ChIP-seq
* Merscope

## Citation

To be added

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
* Filter doublet: `archr_filterdoublet.R`
* Create ArchR project: `archr_createArchRproj.R`
* Peak calling: `archr_callpeaks.R`
* Motif enrichment analysis: `anchr_menr.R`

### 6 Velocity analysis

#### 6.1 RNA velocity 

* Wild-type NC RNA velocity analysis: `scvelo_wt_uncorrected.py`
* foxd3-mutant NC RNA velocity analysis: `scvelo_mut_uncorrected.py`

#### 6.2 MultiVelo 

* Annotate peaks: `01homer_annotpeaks.sh` and `01annot_peaks.R`
* Preprocess ATAC data: `02_1signac_ataccounts.R` and `02_2signac_regionstat.R`
* Link peaks: `03_1run_linkpeaks.sh`
* Preprocess RNA data: `04wtnohox_rna_h5ad.R`
* Prepare ATAC data: `05wtnohox_atac.R`
* Run MultiVelo: `wtnohox_multivelo_uncorrected.py` and `wtnohox_seurat_wnn.R`

#### 6.3 RegVelo

* To be added

### 7 GRN 

#### 7.1 SCENIC+ analysis (`scenicplus`)



#### 7.2 Regulon functional analysis (`regulon_function`)



#### 7.3 GRN dynamics analysis (`GRN_dynamics`)



### 8 ChromBPNet analysis

* Prepare cluster-specific bam: `01prepare_metadata_csv.R`, `02split_bam.sh` and `03merge_bam.sh`
* Preprocess bam: `04prep_bam.sh`
* Preprocess peaks: `04prep_peaks.sh`
* Prepare non peaks and other input: `05chombpnet_prep.sh`
* Run ChromBPNet bias model: `06bias_model.sh`
* Run ChromBPNet bias-facterised model: `07bias_factrised.sh`
* Prepare bigwig: `08pred_bw.sh`
* Compute contribution bigwig: `09contribs_bw.sh`
* TFmodisco de novo motif discovery: `10motif_disc.sh`
* Generate TFmodisco report: `11modisco_report.sh`
* Convert TF modisco output format: `12modisco_convert.sh`
* GIMME cluster motifs: `14gimme_cluster.sh`
* TOMTOM analysis to compare new motifs and existing ones: `14tomtom.sh`
* Footprinting analysis of new motifs: `15footprint_new_motifs.sh`

### 9 AlphaPullDown analysis

* Create individual features: `create_individual_features_SLURM20240202.sh`
* Run multimer: `run_multimer_jobs_SLURM20240202.sh`
* Create notebook: `create_notebook20240202.sh`
* Visualise results: `viz_alphapulldown_rls.Rmd`

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


