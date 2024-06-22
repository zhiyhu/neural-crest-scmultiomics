# Run Sparpros on Wagner et al's organogensis data
# reticulate::use_condaenv("spapros_env")
# Author: Zhiyuan 
# Date: 16 Nov 2023 
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
adata = sc.read_h5ad("ref_mapping/data/Wagner/public_data/WagnerScience2018.processed.h5ad")
adata.var_names = adata.var['_index']
adata.var_names


adata.obs.TimeID.value_counts()
# 24hpf    34750
# 18hpf     6962
# 6hpf      5692
# 10hpf     4280
# 4hpf      4277
# 14hpf     4001
# 8hpf      3568
# Name: TimeID, dtype: int64


import os
os.chdir("STprobe_select/figures/wagner/")

sc.pl.violin(adata, ['tfec','fli1a'], groupby='TimeID', use_raw = False, rotation = 30, save="violinPlots_tfec_fli1a_noLog.png" , )

# Preprocess counts and get highly variable genes
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=1000)

sc.pl.violin(adata, ['tfec', 'fli1a'], groupby='TimeID', use_raw = False, rotation = 30, save="violinPlots_tfec_fli1a.png")

# use case

# list of pre selected genes
preselected_genes = [
    'tfec','tfeb','mitfa','bhlhe40','fli1a','rxraa','rarga','erf','erfl3','elk3','nr2f5','etv2',
    'smarcc1a','ebf3a','nr2f2','ets1','elf1','foxd3','tfap2a'
]

# Get `n_pca_genes` to obtain the top 20 pca genes that are not already in preselected_genes
tmp = sp.se.select_pca_genes(adata, n=adata.n_vars, inplace=False)["selection_ranking"]
n_pca_genes = tmp.loc[~tmp.index.isin(preselected_genes)].sort_values().iloc[:20].max().astype(int)

# create an instance of the ProbesetSelector class
selector = sp.se.ProbesetSelector(
    adata,
    n=70,
    n_pca_genes=n_pca_genes,
    preselected_genes=preselected_genes,
    celltype_key="ClusterName_short",
    verbosity=1,
    n_jobs=-1,
    save_dir="STprobe_select/results/wagner",
)

# select probe set
selector.select_probeset()

with open('STprobe_select/results/wagner/wagner_selector.pkl', 'wb') as outp:
    pickle.dump(selector, outp, pickle.HIGHEST_PROTOCOL)

sp.pl.masked_dotplot(adata, selector, save="masked_dotplot.png")


    
