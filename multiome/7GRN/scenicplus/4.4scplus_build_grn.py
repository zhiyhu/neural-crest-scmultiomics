##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SCENIC+: build GRN
## Zhiyuan Hu
## 28 Nov 2022
## last modified 15 dec 2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus")
# multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/4.4scplus_build_grn.py
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
tmpDir = 'tmp'
r2gdir=work_dir + 'scenicplus/r2g_adj/'
tf2gdir=work_dir + 'scenicplus/tf2g_adj_ss0.2_lr0.01/'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. regions to genes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pickle
infile = open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()

regions_to_genes = {}
for itor in np.arange(0, 17500, 500):
  infile = open(r2gdir + 'r2g_adj' + str(itor) + "_" + str(itor + 500) + '.pkl', 'rb')
  regions_to_genes[str(itor) + "_" + str(itor + 500)] = pickle.load(infile)
  infile.close()

# Adjacency matrix
adj = pandas.concat(regions_to_genes).sort_values(by='importance', ascending=False)

key = 'region_to_gene'
scplus_obj.uns[key] = adj

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. TFs to genes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tfs_to_genes = {}
for itor in np.arange(0, 17600, 200):
  infile = open(tf2gdir + 'tf2g_adj' + str(itor) + "_" + str(itor + 200) + '.pkl', 'rb')
  tfs_to_genes[str(itor) + "_" + str(itor + 200)] = pickle.load(infile)
  infile.close()

# Adjacency matrix
adj = pandas.concat(tfs_to_genes).sort_values(by='importance', ascending=False)

key = 'TF2G_adj'
import time
# log.info('Took {} seconds'.format(time.time() - start_time))
start_time = time.time()
# log.info(f'Adding correlation coefficients to adjacencies.')

ex_matrix = pandas.DataFrame(
    scplus_obj.X_EXP, index=scplus_obj.cell_names, columns=scplus_obj.gene_names)
import scenicplus.TF_to_gene
adj = scenicplus.TF_to_gene._add_correlation(adj, ex_matrix)
adj = scenicplus.TF_to_gene._inject_TF_as_its_own_target(
    TF2G_adj=adj,
    inplace = False,
    ex_mtx = scplus_obj.to_df(layer='EXP'))
# log.info(f'Adding importance x rho scores to adjacencies.')
adj[scenicplus.TF_to_gene.COLUMN_NAME_SCORE_1] = adj[scenicplus.TF_to_gene.COLUMN_NAME_CORRELATION] * \
    adj[scenicplus.TF_to_gene.COLUMN_NAME_WEIGHT]
adj[scenicplus.TF_to_gene.COLUMN_NAME_SCORE_2] = abs(
    adj[scenicplus.TF_to_gene.COLUMN_NAME_CORRELATION]) * abs(adj[scenicplus.TF_to_gene.COLUMN_NAME_WEIGHT])
time.time() - start_time
# log.info('Took {} seconds'.format(time.time() - start_time))
scplus_obj.uns[key] = adj

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Build eGRNs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn

build_grn(scplus_obj,
         min_target_genes = 10,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.03,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         region_to_gene_key = 'region_to_gene',
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = None,
         _temp_dir = os.path.join(tmpDir, 'ray_spill'))


# To access the eGRNs:
import dill
with open(os.path.join(work_dir, 'scenicplus/scplus_obj_grn.pkl'), 'wb') as f:
  dill.dump(scplus_obj, f)

