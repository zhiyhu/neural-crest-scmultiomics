# A Complete Gene Regulatory Network of the Cranial Neural Crest in Zebrafish: single-cell multi-omic and spatial transcriptomic analysis

This repository contains the code for the following analysis, catelogued by the sub-projects:
* [Multiome (RNA+ATAC)](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#10x-multiome-data-analysis)
* [Smart-seq3](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#smart-seq3-data-analysis)
* [Merscope](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#merscope-data-analysis)
* [Perturb-seq](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#perturb-seq-data-analysis)
* [ChIP-seq](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#chip-seq-data-analysis)

## Citation

Zhiyuan Hu, Sarah Mayes, Weixu Wang, José Mariá Santos Perreira, Fabian Theis, Tatjana Sauka-Spengler. A Complete Gene Regulatory Network of the Cranial Neural Crest in Zebrafish (2024).

## 10x multiome data analysis

The analysis scripts below are located in the `multiome` directory.

### 1 Preprocessing and quality control

* Remove ambient RNA by SoupX: `multiome/1preprocessing/soupX_all_samples.R`

* QC and filtering
    * `multiome/1preprocessing/seurat_initialqc_bysample.R`
    <!-- * RNA Doublet filtering: ``
    * Intersect RNA + ATAC good-quality cells: ``
    * Merge objects: `` -->

### 2 Clustering and cell type annotation

* Clustering and annotating all cells
    * Clustering analysis: `multiome/2clustering/1clustering_RNAsoupx_allcells.Rmd`
    * Annotation analysis: `multiome/2clustering/2annotate_RNAsoupx_allcells.Rmd`

* Clustering and annotating NC cells 
    * Clustering analysis: `multiome/2clustering/3clustering_RNAsoupx_NC.Rmd` 
    * Annotation analysis: `multiome/2clustering/4annotate_RNAsoupx_NC.Rmd`

### 3 Genotyping

* Predict genotypes in sample 75% epiboly-4ss: `multiome/3genotyping/amplicon_anlaysis.R`

### 4 Integration

* Use Wagner et al (2018) data to annotate multiome all cells: `multiome/4integration/classification_wagner_multiomeall.R`
* Integrate multiome, Smart-seq3 and Wagner et al (2018) data: `multiome/4integration/cca_NC_wagner2018_multiome_SS3.R`
* Compute inter-group similarity: `multiome/4integration/NC_Wagner2018_CIDER.R`

### 5 ArchR analysis

* Create Arrow files: `multiome/5archr/archr_createArrow.R`
* Filter doublet: `multiome/5archr/archr_filterdoublet.R`
* Create ArchR project: `multiome/5archr/archr_createArchRproj.R`
* Peak calling: `multiome/5archr/archr_callpeaks.R`
* Motif enrichment analysis: `multiome/5archr/anchr_menr.R`

### 6 Velocity analysis

#### 6.1 RNA velocity 

* Wild-type NC RNA velocity analysis: `multiome/6velocity/RNA_velocity/scvelo_wt_uncorrected.py`
* foxd3-mutant NC RNA velocity analysis: `multiome/6velocity/RNA_velocity/scvelo_mut_uncorrected.py`

#### 6.2 MultiVelo 

* Annotate peaks: `multiome/6velocity/multivelo/01homer_annotpeaks.sh` and `multiome/6velocity/multivelo/01annot_peaks.R`
* Preprocess ATAC data: `multiome/6velocity/multivelo/02_1signac_ataccounts.R` and `multiome/6velocity/multivelo/02_2signac_regionstat.R`
* Link peaks: `multiome/6velocity/multivelo/03_1run_linkpeaks.sh`
* Preprocess RNA data: `multiome/6velocity/multivelo/04wtnohox_rna_h5ad.R`
* Prepare ATAC data: `multiome/6velocity/multivelo/05wtnohox_atac.R`
* Run MultiVelo: `multiome/6velocity/multivelo/wtnohox_multivelo_uncorrected.py` and `multiome/6velocity/multivelo/wtnohox_seurat_wnn.R`


### 7 Gene regulatory network (GRN) analysis  

#### 7.1 SCENIC+ analysis (`scenicplus`)

* Reconstruct enhancer-drivern gene regulatory network: scripts in directory `multiome/7GRN/scenicplus`

#### 7.2 Regulon functional analysis (`regulon_function`)

* Centrality analysis: `multiome/7GRN/regulon_function/calculate_centrality_TFonly.py`
* Plot FishEnrichR results for each regulon cluster: `multiome/7GRN/regulon_function/plot_fishEnrichr.R`

#### 7.3 GRN dynamics analysis - SyncReg (`GRN_dynamics_SyncReg`)

* Extract AUC data: `multiome/7GRN/GRN_dynamics_SyncReg/0prepare_auc_mtx.r`
* Extract latent time: `multiome/7GRN/GRN_dynamics_SyncReg/1extract_latent_time.py`
* Preprocess AUC data: `multiome/7GRN/GRN_dynamics_SyncReg/1preprocess_auc_data.R`
* Process regulon data (region-based and gene-based): `multiome/7GRN/GRN_dynamics_SyncReg/2process_regulon_data_rb.R` and `multiome/7GRN/GRN_dynamics_SyncReg/2process_regulon_data_gb.R`
* Compute similarity: `multiome/7GRN/GRN_dynamics_SyncReg/3image_similarity_rb.py` and `multiome/7GRN/GRN_dynamics_SyncReg/3image_similarity_gb.py`
* Perform clustering: `multiome/7GRN/GRN_dynamics_SyncReg/4cluster_image_similarity.R`
* Motif clustering analysis: `multiome/7GRN/GRN_dynamics_SyncReg/tomtom.sh` and `multiome/7GRN/GRN_dynamics_SyncReg/5motif_clustering.R`
* Dendrogram comparison: `multiome/7GRN/GRN_dynamics_SyncReg/6tanglegram_comparison.R`

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

* Quality check, doublet filtering, clustering, finding marker genes: `smartseq3/1preprocessing/qc_8batches.Rmd`

### 2 Clustering

* Subset neural crest cells, clustering, and visualisation: `smartseq3/2clustering/nc_anterior_analysis.Rmd`

### 3 Velocity

* RNA velocity (steady-state model) analysis of Smart-seq3 neural crest data: `smartseq3/3velocity/scvelo_nc01to08_ncwt.ipynb`
*  RNA velocity (dynamic model) analysis: `smartseq3/3velocity/scvelo_ncwt_highRes.ipynb`

## Merscope 

### 1 Gene panel design

* Prepare datasets, run SPARPROS, evaludate the probe set: `merscope/1panel_design`

### 2 Merscope data analysis

* Merge data, clustering: `merscope/2analysis/1merge_run1_run3.ipynb`
* Merscope/snRNA-seq mapping by Tangram: `merscope/2analysis/2run1_tangram_merscopeNC_multiomeNC.ipynb` and `merscope/2analysis/2run3_allCells_mapping.ipynb`
* Spatial visualisation: `merscope/2analysis/3run1_viz_anno20240425.ipynb` and `merscope/2analysis/3run3_viz_anno.ipynb`

## Perturb-seq data analysis

### 1 Preprocessing

* Remove ambient RNA: `perturb-seq/1preprocessing/soupX_all_samples.R`
* Merge data: `perturb-seq/1preprocessing/merge_all_samples.Rmd`
* Doublet identification: `perturb-seq/1preprocessing/doubletfinder_RNAonly.Rmd`

### 2 Clustering and cell type annotation

* Merge SoupX-corrected matrix, cluster all cells and integrate sgRNA information: `perturb-seq/2clustering/1_qc_clustering_allcells_p1top11.Rmd`
* Cluster singlet cells and integrate with multiome data: `perturb-seq/2clustering/2_clustering_allcells_singlet.Rmd`
* Cluster NC cells and integrate with multiome-RNA NC data: `perturb-seq/2clustering/3_clustering_NC_p1top11.Rmd`

### 3 Quantification of in vivo perturbation effects 

* MELD computation and PHATE visualisation: `perturb-seq/3perturbation_effects/1meld_calculateLikelihood.ipynb`
* Impute latent time by RPCA and scVI: `perturb-seq/3perturbation_effects/2impute_by_rpca_results.Rmd` and `perturb-seq/3perturbation_effects/2impute_by_scvi_results.Rmd`
* Quantify perturbation effects: `perturb-seq/3perturbation_effects/3extensive_analysis_meld.Rmd`

### 4 Figure factory

* Effect heatmap, coverage violin plot, PHATE density plot `perturb-seq/4fig_factory/1perturbseq_vis.Rmd`
* Trend plots: `perturb-seq/4fig_factory/2vis_likelihood.Rmd`

## ChIP-seq data analysis

* Snake pipeline: `chip-seq/`
