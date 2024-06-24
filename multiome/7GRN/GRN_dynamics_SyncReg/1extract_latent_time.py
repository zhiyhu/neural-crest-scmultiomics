# Extract latent time
# Zhiyuan Hu
# 29 Jan 2023
# last modified 15 Jan 2024
# reticulate::use_condaenv("multivelo")
# reticulate::repl_python()

import os
import scipy
import numpy as np
import pandas as pd
import numba
pd.__version__
# 1.5.1
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

# import velocyto as vcy

# directories
workdir="multivelo2023dec/"
figdir=workdir+'figures/wt_nohox_uncorrected/'
datadir=workdir+'data/wt_nohox_uncorrected/'
os.chdir(workdir)

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

# run this: multivelo2023dec/code/3multivelo_extract_latentTime.py
adata_result = sc.read_h5ad(datadir+"multivelo_result.h5ad")
mv.velocity_graph(adata_result)
mv.latent_time(adata_result) # adata.uns['iroot'] vkey='velo_u',root_key
adata_result.obs.to_csv("latent_time_grn/data/wtnohox/multivelo_obs.csv")

adata = sc.read_h5ad("velocity/data/scvelo_output/NCwtnohox_scVelo_out_uncorrected.h5ad")
adata.obs.to_csv("latent_time_grn/data/wtnohox/scvelo_obs.csv")
