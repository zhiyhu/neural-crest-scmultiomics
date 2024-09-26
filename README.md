# Single-cell multi-omics, spatial transcriptomics and systematic perturbation decode circuitry of neural crest fate decisions

This repository contains the code for the following analysis, catelogued by the sub-projects:
* [Multiome (RNA+ATAC)](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#10x-multiome-data-analysis)
* [Smart-seq3](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#smart-seq3-data-analysis)
* [Merscope](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#merscope)
* [Perturb-seq](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#perturb-seq-data-analysis)
* [ChIP-seq](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#chip-seq-data-analysis)

## Citation

Zhiyuan Hu*, Sarah Mayes, Weixu Wang, José Mariá Santos Perreira, Fabian Theis, Tatjana Sauka-Spengler*. **Single-cell multi-omics, spatial transcriptomics and systematic perturbation decode circuitry of neural crest fate decisions**. *bioRxiv* (2024). doi: [https://doi.org/10.1101/2024.09.17.613303](https://doi.org/10.1101/2024.09.17.613303)

<img src="https://github.com/zhiyhu/neural-crest-scmultiomics/blob/main/thumb.png" alt="drawing" style="width:600px;"/>

## Resources

* [Neural Crest Cell Browser](https://research.stowers.org/compbio/ucsc_cellbrowser/public_html/): interactively view the 10x multiome (RNA and ATAC) data, Smart-seq3 data, and MERSCOPE data.
* [Interactive Neural Crest GRN](https://research.stowers.org/compbio/NC-GRN_shiny/docs/index.html): interactively inquire the downstream genes and regions of individual regulons or interactions between two regulons.

## 10x Multiome data analysis

The analysis scripts below are located in the `multiome` directory.

### 1 Preprocessing and quality control

Directory: `multiome/1preprocessing/`

* Remove ambient RNA by [SoupX](https://github.com/constantAmateur/SoupX): `0soupX_all_samples.R`

* QC and filtering
    * RNA data QC: `1seurat_initialqc_bysample.R`
    * RNA doublet filtering by [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder): `2doubletfinder_RNA/*.Rmd`
    * ATAC data QC and doublet filtering: see [ArchR analysis](https://github.com/zhiyhu/neural-crest-scmultiomics?tab=readme-ov-file#5-archr-analysis)
    * Intersect RNA + ATAC good-quality cells: `4intersect_rna_atac_singlet.R`
    * Merge objects: `5merge_seuobj_rna.R`

### 2 Clustering and cell type annotation 

Directory: `multiome/2clustering/`

* Clustering and annotating all cells
    * Clustering analysis: `1clustering_RNAsoupx_allcells.Rmd`
    * Annotation analysis: `2annotate_RNAsoupx_allcells.Rmd`

* Clustering and annotating NC cells 
    * Clustering analysis: `3clustering_RNAsoupx_NC.Rmd` 
    * Annotation analysis: `4annotate_RNAsoupx_NC.Rmd`

### 3 Genotyping

* Predict genotypes in sample 75% epiboly-4ss: `multiome/3genotyping/amplicon_anlaysis.R`

### 4 Integration

Directory: `multiome/4integration/`

* Use [Wagner et al (2018)](https://www.science.org/doi/10.1126/science.aar4362) data to annotate multiome all cells: `classification_wagner_multiomeall.R`
* Integrate multiome, Smart-seq3 and Wagner et al (2018) data: `cca_NC_wagner2018_multiome_SS3.R`
* Compute inter-group similarity by [CIDER](https://github.com/zhiyhu/CIDER): `NC_Wagner2018_CIDER.R`

### 5 ArchR analysis

Directory: `multiome/5archr/`

* Create Arrow files by [ArchR](https://www.archrproject.com): `archr_createArrow.R`
* Filter doublet: `archr_filterdoublet.R`
* Create ArchR project: `archr_createArchRproj.R`
* Peak calling: `archr_callpeaks.R`
* Motif enrichment analysis: `anchr_menr.R`

### 6 Velocity analysis

#### 6.1 RNA velocity (`multiome/6velocity/RNA_velocity/`)

* Wild-type NC RNA velocity analysis: `scvelo_wt_uncorrected.py`
* foxd3-mutant NC RNA velocity analysis: `scvelo_mut_uncorrected.py`

#### 6.2 MultiVelo (`multiome/6velocity/multivelo/`)

* Annotate peaks: `01homer_annotpeaks.sh` and `01annot_peaks.R`
* Preprocess ATAC data: `02_1signac_ataccounts.R` and `02_2signac_regionstat.R`
* Link peaks: `03_1run_linkpeaks.sh`
* Preprocess RNA data: `04wtnohox_rna_h5ad.R`
* Prepare ATAC data: `05wtnohox_atac.R`
* Run [MultiVelo](https://multivelo.readthedocs.io): `wtnohox_multivelo_uncorrected.py` and `wtnohox_seurat_wnn.R`


### 7 Gene regulatory network (GRN) analysis  

#### 7.1 SCENIC+ analysis

* Reconstruct enhancer-drivern gene regulatory network by [SCENIC+](https://scenicplus.readthedocs.io): scripts in directory `multiome/7GRN/scenicplus`

#### 7.2 Regulon functional analysis

Directory: `multiome/7GRN/regulon_function`

* Regulon clustering analysis: `cluster_regulons_regionBased.py`
* Centrality analysis: `calculate_centrality_TFonly.py`
* Visualisation of clustering and centrality analysis results: `regulon_clustering_centrality.Rmd`
* Visualisation of [FishEnrichR](https://maayanlab.cloud/FishEnrichr/) results: `plot_fishEnrichr.R`

#### 7.3 GRN dynamics analysis - SyncReg 

Directory: `multiome/7GRN/GRN_dynamics_SyncReg`

* Extract AUC data: `0prepare_auc_mtx.r`
* Extract latent time: `1extract_latent_time.py`
* Preprocess AUC data: `1preprocess_auc_data.R`
* Process regulon data (region-based and gene-based): `2process_regulon_data_rb.R` and `2process_regulon_data_gb.R`
* Compute similarity: `3image_similarity_rb.py` and `3image_similarity_gb.py`
* Perform clustering: `4cluster_image_similarity.R`
* Motif clustering analysis: `tomtom.sh` and `5motif_clustering.R`
* Dendrogram comparison: `6tanglegram_comparison.R`

### 8 ChromBPNet analysis

* Prepare cluster-specific bam: `01prepare_metadata_csv.R`, `02split_bam.sh` and `03merge_bam.sh`
* Preprocess bam: `04prep_bam.sh`
* Preprocess peaks: `04prep_peaks.sh`
* Prepare non peaks and other input: `05chombpnet_prep.sh`
* Run [ChromBPNet](https://github.com/kundajelab/chrombpnet) bias model: `06bias_model.sh`
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
* Run multimer via [AlphaPullDown](https://github.com/KosinskiLab/AlphaPulldown): `run_multimer_jobs_SLURM20240202.sh`
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

* Prepare datasets, run [SPAPROS](https://spapros.readthedocs.io), evaluate the probe set: `merscope/1panel_design`

### 2 Merscope data analysis

Directory: `merscope/2analysis`

* Merge data, clustering: `1merge_run1_run3.ipynb`
* Merscope/snRNA-seq mapping by [Tangram](https://tangram-sc.readthedocs.io): `2run1_tangram_merscopeNC_multiomeNC.ipynb` and `2run3_allCells_mapping.ipynb`
* Spatial visualisation: `3run1_viz_anno20240425.ipynb` and `3run3_viz_anno.ipynb`

## Perturb-seq data analysis

### 1 Preprocessing

Directory: `perturb-seq/1preprocessing/`

* Remove ambient RNA: `soupX_all_samples.R`
* Merge data: `merge_all_samples.Rmd`
* Doublet identification: `doubletfinder_RNAonly.Rmd`

### 2 Clustering and cell type annotation

Directory: `perturb-seq/2clustering/`

* Merge SoupX-corrected matrix, cluster all cells and integrate sgRNA information: `1_qc_clustering_allcells_p1top11.Rmd`
* Cluster singlet cells and integrate with multiome data: `2_clustering_allcells_singlet.Rmd`
* Cluster NC cells and integrate with multiome-RNA NC data: `3_clustering_NC_p1top11.Rmd`

### 3 Quantification of in vivo perturbation effects 

Directory: `perturb-seq/3perturbation_effects/`

* [MELD](https://github.com/KrishnaswamyLab/MELD) computation and PHATE visualisation: `1meld_calculateLikelihood.ipynb`
* Impute latent time by [RPCA](https://satijalab.org/seurat/articles/integration_rpca.html) and [scVI](https://docs.scvi-tools.org/en/stable/tutorials/index.html): `2impute_by_rpca_results.Rmd` and `2impute_by_scvi_results.Rmd`
* Quantify perturbation effects: `3extensive_analysis_meld.Rmd`

### 4 Figure factory 

Directory: `perturb-seq/4fig_factory/`

* Effect heatmap, coverage violin plot, [PHATE](https://phate.readthedocs.io/en/stable/tutorial.html) density plot `1perturbseq_vis.Rmd`
* Trend plots: `2vis_likelihood.Rmd`

## ChIP-seq data analysis

* Snake pipeline: `chip-seq/`
