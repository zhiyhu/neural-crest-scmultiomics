# Run scVelo on anterior NC cells; MUT; uncorrected
# reticulate::use_condaenv("cellrank_new")
# reticulate::repl_python()
# Zhiyuan Hu
# Created 23 Feb 2023
# Last modified 30 March 2023

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
# import seaborn as sns
import matplotlib as plt
import os
import loompy

scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
os.chdir('multiome/analysis_newref/velocity/')

proj="scvelo_mut_uncorrected"

adata_seurat = scv.read("data/seurat/seu_RNAsoupx_NC.h5ad")
adata = scv.read("data/rds/NC_RNAvelocyto_uncorrected.h5ad", cache=True)
adata = scv.utils.merge(adata, adata_seurat)
adata

# # frequency of cluster+genotype
# pd.crosstab(adata.obs['seurat_clusters'],adata.obs['genotype_new']).stack().reset_index(name='Freq')

adata.obs["genotype_new"].value_counts()

####### subset data
adata = adata[adata.obs['seurat_clusters'].isin([4,7,9,12,13,14,15,16,17,18,19,20]),:].copy()
adata = adata[adata.obs['genotype_new'].isin(['dp']),:].copy()

adata.obs['sample'] = adata.obs['sample'].astype("category")
scv.pl.proportions(adata, groupby = 'sample', save="../velocity/figures/"+proj+"/scvelo_proportions_sample.pdf")

### Scanpy QC
adata.var_names=adata.var["features"]
sc.pl.highest_expr_genes(adata, n_top=20)

###  Preprocessing the data
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata,enforce=True)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
adata

# compute neighbours
scv.pp.neighbors(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=10)

### Estimate RNA velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

### Project the velocities
scv.pl.velocity_embedding(adata, basis='tsne', color='seurat_clusters', arrow_length=3, arrow_size=2, dpi = 120,
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding.pdf")
scv.pl.velocity_embedding_grid(adata, basis='tsne', color='seurat_clusters',arrow_length=3, 
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding_grid.pdf")
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='seurat_clusters', density = 5, dpi = 120, 
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding_stream.png")
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='stage', density = 5, dpi = 120, 
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding_stream_stage.png")

#-------------------------------#
# plot using new colour palette #
#-------------------------------#
colours_df=pd.read_csv("/Filers/home/z/zhu/t1data/multiome/analysis_newref/clustering/figures/for_pre/df_coloursUserd.csv")
# Convert the first and third columns to a dictionary
col_dict = dict(zip(colours_df['cell_type'].tolist(), colours_df['colour'].tolist()))
# velocity embedding arrows
scv.pl.velocity_embedding(adata, basis='tsne', color='cell_type', 
arrow_length=3, arrow_size=2, dpi = 120, palette=col_dict,
legend_loc='right',
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding_newPalette.pdf")

# velocity embedding streams
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='cell_type', 
density = 5,  dpi = 120, palette=col_dict,
legend_loc='right',
save="../velocity/figures/"+proj+"/scvelo_velocity_embedding_stream_newPalette.pdf")

### Speed and coherence

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, basis = "tsne",cmap='coolwarm', perc=[5, 95], 
save="../velocity/figures/"+proj+"/scvelo_scatter_speed_confidence.pdf")

df = adata.obs.groupby('seurat_clusters')[keys].mean().T
df

df.style.background_gradient(cmap='coolwarm', axis=1).to_excel("../velocity/figures/"+proj+"/scvelo_speed_confidence.xlsx", engine='openpyxl')

### Velocities and pseudotime

scv.pl.velocity_graph(adata, basis='tsne', threshold=.3,dpi = 120, 
save="../velocity/figures/"+proj+"/scvelo_velocity_graph_threshold.1.pdf")

# the idx of starting cell

starting_idx=np.where([adata.obs_names == 's8_TTGTCCCAGGTATTGC-1'])[1][0]
x, y = scv.utils.get_cell_transitions(adata, basis='tsne',starting_cell=starting_idx)
ax = scv.pl.velocity_graph(adata, c='lightgrey', basis='tsne', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, basis='tsne', c='ascending', 
cmap='gnuplot', ax=ax)
plt.pyplot.savefig("../velocity/figures/"+proj+"/scvelo_cell_transitions_gnuplot.pdf")  


# scv.tl.terminal_states(adata)
# scv.pl.scatter(adata, color=['root_cells', 'end_points'], basis='tsne')

scv.tl.velocity_pseudotime(adata, n_dcs=24, root_key=starting_idx)
scv.pl.scatter(adata, basis='tsne', color='velocity_pseudotime', cmap='gnuplot',
save="../velocity/figures/"+proj+"/scvelo_velocity_pseudotime.pdf")

### PAGA velocity graph
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='seurat_clusters',root_key=starting_idx,
            minimum_spanning_tree = False, threshold_root_end_prior = 0.5)
    
df = scv.get_df(adata, 'paga/transitions_confidence').T
df.style.background_gradient(cmap='Blues').format('{:.2g}').to_excel("../velocity/figures/"+proj+"/scvelo_paga_transitions_confidence.xlsx", engine='openpyxl')


scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_paga.pdf")

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.04,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_paga_threshold0.04.pdf")

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.03,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_paga_threshold0.03.pdf")

## new palette
scv.tl.paga(adata, groups='cell_type',root_key=starting_idx,
            minimum_spanning_tree = False, threshold_root_end_prior = 0.5)
    
scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0, palette = col_dict,
            min_edge_width=1, node_size_scale=0.5, 
            save="../velocity/figures/"+proj+"/scvelo_paga_newPalatte.pdf")

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.04, palette = col_dict,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_paga_threshold0.04_newPalatte.pdf")
 
scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.03, palette = col_dict,
            min_edge_width=1, node_size_scale=0.5, legend_loc='right',
            save="../velocity/figures/"+proj+"/scvelo_paga_threshold0.03_newPalatte.pdf")


sc.pl.tsne(adata, color= ["elavl3"], use_raw = False, save = "/../"+proj+"/elval3_expr_tsne.pdf")


###-------------------------------
### Save data

## Fix a bug of anndata: https://github.com/theislab/scvelo/issues/255
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
del(adata.var['_index'])
adata.write_h5ad("data/scvelo_output/"+proj+"scVelo_out.h5ad")

###-------------------------------
## Interpret the velocity

adata.var_names=adata.var["features"]
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
for x in ['4', '7', '9', '12', '13', '14', '15', '16', '17', '18', '19', '20']:
  scv.pl.velocity(adata, df[x][0:6], basis = "tsne", color="seurat_clusters", ncols=2,
  save="../velocity/figures/"+proj+"/NC_velocity_Cluster"+x+"_top6.pdf")

###--------------------------------
### dynamical modelling
###--------------------------------

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

## add the new colour panel here
colours_df=pd.read_csv("/Filers/home/z/zhu/t1data/multiome/analysis_newref/clustering/figures/for_pre/df_coloursUserd.csv")
# Convert the first and third columns to a dictionary
col_dict = dict(zip(colours_df['cell_type'].tolist(), colours_df['colour'].tolist()))

# embedding stream with labels
scv.pl.velocity_embedding_stream(adata, basis='tsne',color = "cell_type", 
density = 5, dpi = 120, palette=col_dict,
save="../velocity/figures/"+proj+"/scvelo_dynamical_velocity_embedding_stream.pdf")

# embedding stream with no label over the figure
scv.pl.velocity_embedding_stream(adata, basis='tsne',color = "cell_type", 
density = 5, dpi = 120, palette=col_dict,legend_loc='right',
save="../velocity/figures/"+proj+"/scvelo_dynamical_velocity_embedding_stream_noLabel.pdf")

### PAGA velocity graph
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='cell_type',root_key=starting_idx,
            minimum_spanning_tree = False, threshold_root_end_prior = 0.5)
    
df = scv.get_df(adata, 'paga/transitions_confidence').T
df.style.background_gradient(cmap='Blues').format('{:.2g}').to_excel('../velocity/figures/"+proj+"/scvelo_dynamical_paga_transitions_confidence.xlsx', engine='openpyxl')

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_dynamical_paga.pdf")

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.04,
            min_edge_width=1, node_size_scale=0.5,
            save="../velocity/figures/"+proj+"/scvelo_dynamical_paga_threshold0.04.pdf")

scv.pl.paga(adata, basis='tsne', size=50, alpha=.1,threshold=0.03,
            min_edge_width=1, node_size_scale=0.5, legend_loc='right',
            save="../velocity/figures/"+proj+"/scvelo_dynamical_paga_threshold0.03.pdf")

# latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, basis="tsne",
save="../velocity/figures/"+proj+"/scvelo_dynamical_latent_time.pdf")
