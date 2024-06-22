# Run Sparpros on SS3 NC data
# reticulate::use_condaenv("spapros_env")
# Author: Zhiyuan 
# Date: 15 Nov 2023 
# Last modified: 17 nov 2023

import sys
print(sys.path)

import numpy as np
import scanpy as sc
import spapros as sp

sc.settings.verbosity = 0
sc.logging.print_header()
print(f"spapros=={sp.__version__}")

# load and process data

adata = sc.read("integration/data/ss3_ncall.h5ad")
# adata.var_names = adata.var['_index']
adata.var_names

adata

import os
os.chdir("STprobe_select/figures/ss3/ncall/")
sc.pl.violin(adata, ['tfec', 'fli1a'], groupby='cell_type', use_raw = False, rotation = 30, save="violinPlots_tfec_fli1a_noLog.png" )

# Preprocess counts and get highly variable genes
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1000)

sc.pp.filter_genes(adata, min_cells=50)

adata

import os
os.chdir("STprobe_select/figures/ss3/ncall/")
sc.pl.violin(adata, ['tfec', 'fli1a'], groupby='cell_type', use_raw = False, rotation = 30, save="violinPlots_tfec_fli1a_log1p.png")

# use case

# list of pre selected genes
preselected_genes = [
    'tfec','tfeb','mitfa','bhlhe40','fli1a','rxraa','rarga','erf','erfl3','elk3','nr2f5','etv2',
    'smarcc1a','ebf3a','nr2f2','ets1','elf1','foxd3','tfap2a',
    'foxc1a','postnb','dlx2a','grem2b','itga3a','cadm1a','erbb3b','hoxc3a'
]

# Get `n_pca_genes` to obtain the top 20 pca genes that are not already in preselected_genes
tmp = sp.se.select_pca_genes(adata, n=adata.n_vars, inplace=False)["selection_ranking"]
n_pca_genes = tmp.loc[~tmp.index.isin(preselected_genes)].sort_values().iloc[:20].max().astype(int)
n_pca_genes
# create an instance of the ProbesetSelector class
selector = sp.se.ProbesetSelector(
    adata,
    n=70,
    n_pca_genes=n_pca_genes,
    preselected_genes=preselected_genes,
    celltype_key="cell_type",
    verbosity=1,
    n_jobs=-1,
    save_dir="STprobe_select/results/ss3"
)

# select probe set
selector.select_probeset()

first_selection = selector.probeset[selector.probeset["selection"]].index.tolist()

selector.probeset.index[selector.probeset.selection]

sp.pl.masked_dotplot(adata, selector, ct_key = "cell_type",save="masked_dotplot.png")

import dill as pickle
with open('STprobe_select/results/ss3/ss3_ncall_selector.pkl', 'wb') as outp:
   pickle.dump(selector, outp, pickle.HIGHEST_PROTOCOL)

    
