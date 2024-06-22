##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Scanpy process for SCENIC+:  All NC clusters and NO downsampling
## Zhiyuan Hu
## 22 Nov 2022
## last modified 21 Jan 2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reticulate::use_condaenv("scenicplus")

# n_cells = 5000 # number of cells to keep
work_dir = 'multiome/analysis_newref/GRN_scenicplus/ncall/'
tmpDir = 'tmp'

import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt

#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))
os.chdir(os.path.join(work_dir, 'scRNA'))

###### read data  ###### 
adata = sc.read_h5ad('multiome/analysis_newref/velocity/data/seurat/seu_RNAsoupx_NC.h5ad')
# fix raw issue
adata_tmp = adata.raw.to_adata()
adata_tmp.var_names = adata.var['features']
adata.raw = adata_tmp
# add var_names
adata.var_names = adata.var['features']
adata

######## QC  ###### 
sc.pl.scatter(adata, x='nCount_RNA', y='percent.mt', show=True, save =  "adata_countVSmt.pdf")
sc.pl.scatter(adata, x='nCount_RNA', y='nFeature_RNA', show = True, save ="adata_countVSfeature.pdf")

###### tSNE ###### 
import os
os.chdir(os.path.join(work_dir, 'scRNA'))
sc.pl.tsne(adata, color = 'cell_type', save =  "adata_cell_annot.pdf")

#### Modifiy foxd3 expression ###### 
adata.obs['factor'] = 1
adata.obs['factor'][adata.obs['genotype_new'] == "dp"] = 10

adata.var_names = adata.var['_index']
idx=np.where(adata.var_names == 'foxd3')
idx
adata.X[:,int(idx[0])] = (adata.X.tocoo().tocsc()[:,int(idx[0])].todense().ravel()/adata.obs['factor'].ravel())
sc.pl.tsne(adata, color=['foxd3-mCherry','foxd3','genotype_new'], use_raw=False, save =  "adata_modified_foxd3expr_beforelog.pdf")

# save unprocessed data
expr_mat = adata.X
expr_mat.index = adata.obs_names

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pl.tsne(adata, color=['foxd3-mCherry','foxd3','genotype_new'], use_raw=False, save =  "adata_modified_foxd3expr.pdf")
sc.pl.tsne(adata, color=['foxd3-mCherry','foxd3','genotype_new'], use_raw=True, save =  "adata_modified_foxd3expr_raw.pdf")

###### Create new anndata ###### 

import anndata
rna_anndata = anndata.AnnData(X=expr_mat, 
  dtype= "float64",
  obs = adata.obs,
  var = adata.var,
  obsm = adata.obsm,
  varm = adata.varm)
rna_anndata.obs = adata.obs

###### save data ###### 
# Save
del rna_anndata.var['_index']
rna_anndata.write_h5ad(os.path.join(work_dir,"scRNA/adata_ncall.h5ad"))
adata.var.to_csv(os.path.join(work_dir,"scRNA/adata_rna_genenames.csv"))
