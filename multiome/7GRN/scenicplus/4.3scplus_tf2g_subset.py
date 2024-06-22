##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inferring eGRNs using SCENIC+ 
## Zhiyuan Hu
## 26 Nov 2022
## last modified 15 dec 2023
##
## Env: scenicplus
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus")
# Refer to https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html

import dill
import scanpy as sc
import numpy as np
import os
import anndata
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/'
tmpDir = 'tmp3/'

idx_start = sys.argv[1]
idx_end = sys.argv[2]

idx_start = int(idx_start)
idx_end = int(idx_end)

import pickle
infile = open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()

# avoid the last one over scplus_obj max gene numbers
idx_end = min(idx_end, len(scplus_obj.gene_names))

print('Index start: ', idx_start)
print('Index end: ', idx_end)
genes_subset = scplus_obj.gene_names.tolist()[idx_start:idx_end]
print('Gene subset:', genes_subset)

# from scenicplus.TF_to_gene import *
tf_file = 'multiome/analysis_newref/GRN_scenicplus/data/motif/dr11_AnimalTFDB4_DanioCode_TF_lst_curated20230122.txt'

SGBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 5000,  # can be arbitrarily large
    'max_features': 0.1,
    'subsample': 0.2
}
method_params = [
  'GBM',      # regressor_type
  SGBM_KWARGS  # regressor_kwargs
]

print(method_params)
sys.path.append("multiome/analysis_newref/GRN_scenicplus/code/custom_functions")
import functions_tf2g_subset
tf2g_adj = functions_tf2g_subset.calculate_TFs_to_genes_relationships_my(scplus_obj,
                    tf_file = tf_file,
                    method = 'GBM',
                    method_params = method_params,
                    genes = genes_subset,
                    key= 'TF2G_adj')

# Save
if not os.path.exists(os.path.join(work_dir, 'scenicplus/tf2g_adj_ss0.2_lr0.01')):
    os.makedirs(os.path.join(work_dir, 'scenicplus/tf2g_adj_ss0.2_lr0.01'))

# Save
import pickle
with open(work_dir + 'scenicplus/tf2g_adj_ss0.2_lr0.01/tf2g_adj' + str(idx_start) + '_' + str(idx_end) + '.pkl', 'wb') as f:
  pickle.dump(tf2g_adj, f)


