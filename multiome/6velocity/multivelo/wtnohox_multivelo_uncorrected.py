# Run multivelo on the data anterior NC-mutant
# /ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/code/wtnohox_multivelo_uncorrected.py
# Zhiyuan Hu
# 23 Feb 2023
# last modified 29 dec 2023
# reticulate::use_condaenv("multivelo_new")
# reticulate::repl_python()

import os
import scipy
import numpy as np
import pandas as pd
import numba
pd.__version__
# 1.5.2
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

# import velocyto as vcy

# directories
workdir="/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/"
figdir=workdir+'figures/wt_nohox_uncorrected/'
datadir=workdir+'data/wt_nohox_uncorrected/'
os.chdir(workdir)

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading in unspliced and spliced counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read the spliced/unspliced reads
adata_rna = scv.read(workdir + "/data/wt_nohox/rna_adata_uncorrected.h5ad")
adata_rna
# AnnData object with n_obs × n_vars = 6970 × 24773

# Load cell barcodes
cell_barcodes = pd.read_csv(workdir + '/data/wt_nohox/cell_barcodes.txt', sep='\t', header=None)
adata_rna = adata_rna[cell_barcodes[0],].copy()
adata_rna_raw = adata_rna

sc.pp.filter_cells(adata_rna, min_counts=1000)
sc.pp.filter_cells(adata_rna, max_counts=20000)
adata_rna
# AnnData object with n_obs × n_vars = 6957 × 24773

scv.pp.filter_genes(adata_rna, min_shared_counts=10)
scv.pp.normalize_per_cell(adata_rna, enforce=True)
scv.pp.filter_genes_dispersion(adata_rna, n_top_genes=1000)
scv.pp.log1p(adata_rna)
adata_rna
# AnnData object with n_obs × n_vars = 6957 × 1000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preprocessing the ATAC counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
adata_atac = sc.read("/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox/atac_adata_cistopicPeaks.h5ad")
adata_atac
# AnnData object with n_obs × n_vars = 6970 × 275579

# We aggregate peaks around each gene as well as those that have high correlations with promoter peak or gene expression.
# Peak annotation contains the metadata for all peaks.
# Feature linkage contains pairs of correlated genomic features.
import sys
sys.path.append("/ceph/project/tsslab/zhu/multiome/R/multivelo/code/multivelo_modified")
import aggregate_peaks_10x

adata_atac = aggregate_peaks_10x.aggregate_peaks_10x_my(adata_atac,
                                    peak_annot_file = '/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/01peaks/multivelo_peakanno_scplus.txt',
                                    linkage_file = '/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/03feature_linkage/signac_feature_linkage_all.bedpe',
                                    verbose=True)
adata_atac
# View of AnnData object with n_obs × n_vars = 6970 × 15436

# Let's examine the total count distribution and remove outliers.
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000))     
plt.savefig(figdir+ 'atac_total_counts_hist.png')

max(adata_atac.X.sum(1))
min(adata_atac.X.sum(1))

sc.pp.filter_cells(adata_atac, min_counts=1000)
sc.pp.filter_cells(adata_atac, max_counts=60000)

# We normalize aggregated peaks with TF-IDF.
mv.tfidf_norm(adata_atac)

adata_atac
# AnnData object with n_obs × n_vars = 5698 × 15436
#     obs: 'n_counts'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Finding shared barcodes and features between RNA and ATAC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var._index, adata_atac.var_names))
len(shared_cells), len(shared_genes)
# (5692, 853)

adata_rna.var_names = adata_rna.var._index

adata_rna = adata_rna_raw.copy()
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

# renormalise RNA data again after filtering???
scv.pp.normalize_per_cell(adata_rna, enforce=True)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=30)

adata_rna.obs['cell_type'] = adata_rna.obs['cell_type'].astype('category')

# Reorder the categories for color consistency with the manuscript.
all_clusters = ['NPB_nohox', 'NPB_nohox_cycling',
                'dNC_nohox','dNC_nohox_cycling',
                'mNC_nohox', 'mNC_arch1', 'mNC_head_mesenchymal',
                'Pigment_sox6_high','Pigment_gch2_high']
adata_rna.obs['cell_type'] = adata_rna.obs['cell_type'].cat.reorder_categories(all_clusters)
adata_rna.obs['stage'] = adata_rna.obs['stage'].astype('category')
adata_rna.obs['stage'] = adata_rna.obs['stage'].cat.reorder_categories(['epiboly-4ss','6-10ss','12-16ss','18-22ss'])

# umap
os.chdir(figdir)
scv.tl.umap(adata_rna)
scv.pl.umap(adata_rna, color = "cell_type", save = 'multivelo_wtnohox_corrected.pdf')

# --------------------------------------------- #
# Smoothing gene aggregagted peaks by neighbors #
# --------------------------------------------- #

# Write out filtered cells and prepare to run Seurat WNN --> R script can be found on Github.
os.chdir(datadir)
adata_rna.obs_names.to_frame().to_csv('filtered_cells.txt', header=False, index=False)

# run /ceph/home/z/zhu/t1data/multiome/analysis_newref/multivelo2023dec/code/wtnohox_seurat_wnn.R

# Read in Seurat WNN neighbors.
os.chdir('/ceph/project/tsslab/zhu/multiome/analysis_newref/multivelo2023dec/data/wt_nohox_uncorrected/')
nn_idx = np.loadtxt("seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("seurat_wnn/nn_dist.txt", delimiter=',')
nn_cells = pd.Index(pd.read_csv("seurat_wnn/nn_cells.txt", header=None)[0])
os.chdir(datadir)

# # Make sure cell names match.
np.all(nn_cells == adata_atac.obs_names)
mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

adata_atac
# AnnData object with n_obs × n_vars =  5692 × 853

# ---------------------------------- #
# Running multi-omic dynamical model #
# ---------------------------------- #

# This will take a while. Parallelization is high recommended.
adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac, 
                                         max_iter=5, 
                                         init_mode="invert", 
                                         verbose=True, 
                                         parallel=True, 
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500, 
                                         extra_color_key='cell_type'
                                        )


adata_result.__dict__['_raw'].__dict__['_var'] = adata_result.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del adata_result.var["_index"]
adata_result.write(datadir+"multivelo_result.h5ad")

# adata_result = sc.read_h5ad(datadir+"multivelo_result.h5ad")

# plot
mv.pie_summary(adata_result)
plt.savefig(figdir+'multivelo_pie_summary.pdf')

# plot
mv.switch_time_summary(adata_result)
plt.savefig(figdir+'multivelo_switch_time_summary.pdf')

# plot
mv.likelihood_plot(adata_result)
plt.savefig(figdir+'multivelo_likelihood_plot.pdf')

# ----------------------------------------- #
# Computing velocity stream and latent time #
# ----------------------------------------- #

mv.velocity_graph(adata_result)
mv.latent_time(adata_result) # adata.uns['iroot'] vkey='velo_u',root_key

scv.tl.velocity_pseudotime(adata_result,vkey='velo_s'+'_norm')

import seaborn as sns
fig, ax = plt.subplots(figsize=(4, 5))
ax = sns.boxplot(data=adata_result.obs, x="stage", y="latent_time", ax =ax)
fig = ax.get_figure()
fig.savefig(figdir+"stage_vs_latenttime.pdf")
    
import seaborn as sns
fig, ax = plt.subplots(figsize=(4, 5))
ax = sns.boxplot(data=adata_result.obs, x="stage", y="velo_s_norm_pseudotime", ax =ax)
fig = ax.get_figure()
fig.savefig(figdir+"stage_vs_pseudotime.pdf")
 
mv.velocity_embedding_stream(adata_result, basis='umap', color='cell_type')
plt.savefig(figdir+'multivelo_velocity_embedding_stream.png')

# a new figure for the presentation
## add the new colour panel here
colours_df=pd.read_csv("/Filers/home/z/zhu/t1data/multiome/analysis_newref/clustering/figures/for_pre/df_coloursUserd.csv")
# Convert the first and third columns to a dictionary
col_dict = dict(zip(colours_df['cell_type'].tolist(), colours_df['colour'].tolist()))

mv.velocity_embedding_stream(adata_result, basis='umap', color='cell_type',legend_loc='right',
density = 3, dpi = 120, palette=col_dict, save = figdir+'multivelo_velocity_embedding_stream_newColourDict.png')
# plt.savefig(figdir+'multivelo_velocity_embedding_stream_newColourDict.png')

# scatter latent time
scv.pl.scatter(adata_result, color='latent_time', basis='umap', color_map='gnuplot', size=80)
plt.savefig(figdir+'multivelo_latent_time_scatter.pdf')

scv.pl.scatter(adata_result, color='velo_s_norm_pseudotime', basis='umap', color_map='gnuplot', size=80)

adata_result.var.to_csv(datadir+"multivelo_var_data.csv")

# ----------------------------------------- #
# Let’s examine some example genes.         #
# ----------------------------------------- #
gene_list = ['col11a1a','col1a1a','ednraa','elavl3','fgfr4','hoxb3a',
'sox5','sox6','sox13','snai1b','sox9b','tfap2b','tfap2e','zeb2b','bhlhe40']

# We can plot accessbility and expression against gene time.
# Accessibility/expression by gene time, colored by the four potential states.
# The solid black curve indicates anchors.
mv.dynamic_plot(adata_result, gene_list, color_by='state', axis_on=False, frame_on=False)
plt.savefig(figdir+'multivelo_dynamic_plot.pdf')

# We can plot velocity against gene time
# Velocity by gene time, colored by the four potential states.
# The solid black curve indicates anchors.
mv.dynamic_plot(adata_result, gene_list, color_by='state', by='velocity', axis_on=False, frame_on=False)
plt.savefig(figdir+'multivelo_dynamic_plot_colour_by_state.pdf')

# We can examine accessibility and expression against globally shared latent time.
# Accessibility/expression by global latent time, colored by cell type assignments.
# The solid black curve indicates the mean.
mv.dynamic_plot(adata_result, gene_list, color_by='cell_type', gene_time=False, axis_on=False, frame_on=False)
plt.savefig(figdir+'multivelo_dynamic_plot_colour_by_celltype.pdf')

# Unspliced-spliced phase portraits, colored by celltype.
mv.scatter_plot(adata_result, gene_list, color_by='cell_type', by='us', axis_on=False, frame_on=False)
plt.savefig(figdir+'multivelo_plot_us_colour_by_celltype.pdf')

# Chromatin-unspliced phase portraits, colored by celltype.
mv.scatter_plot(adata_result, gene_list, color_by='cell_type', by='cu', axis_on=False, frame_on=False)
plt.savefig(figdir+'multivelo_plot_cu_colour_by_celltype.pdf')
