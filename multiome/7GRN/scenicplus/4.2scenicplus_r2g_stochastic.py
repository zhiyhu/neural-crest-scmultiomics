##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inferring eGRNs using SCENIC+: all NC cells
## Zhiyuan Hu
## 21 Nov 2022
## last modified 15 Dec 2023
##
## Env: scenicplus
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus")
# Refer to https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html

import dill
import scanpy as sc
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
work_dir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/'
tmpDir = '/ceph/home/z/zhu/t1data/tmp3'
ray_n_cpu = 48

### start and end index
idx_start = sys.argv[1]
idx_end = sys.argv[2]

idx_start = int(idx_start)
idx_end = int(idx_end)

import pickle
infile = open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()
# idx_start=0
# idx_end=len(scplus_obj.gene_names)

# avoid the last one over scplus_obj max gene numbers
idx_end = min(idx_end, len(scplus_obj.gene_names))

print('Index start: ', idx_start)
print('Index end: ', idx_end)

genes_subset = scplus_obj.gene_names.tolist()[idx_start:idx_end]
print('Gene subset:', genes_subset)

# from scenicplus.enhancer_to_gene import calculate_regions_to_genes_relationships
sys.path.append('/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/code/custom_functions')
import enhancer_to_gene_my 

SGBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 5000,  # can be arbitrarily large
    'max_features': 0.1,
    'subsample': 0.5
}
print(SGBM_KWARGS)

print("Calculate_regions_to_genes_relationships starts")
r2g_adj=enhancer_to_gene_my.calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = None,
                    genes= genes_subset,
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = SGBM_KWARGS, inplace = False)

print("Calculate_regions_to_genes_relationships done")


# Save
if not os.path.exists(os.path.join(work_dir, 'scenicplus/r2g_adj')):
    os.makedirs(os.path.join(work_dir, 'scenicplus/r2g_adj'))

import pickle
with open(work_dir + 'scenicplus/r2g_adj/r2g_adj' + str(idx_start) + '_' + str(idx_end) + '.pkl', 'wb') as f:
  pickle.dump(r2g_adj, f)
