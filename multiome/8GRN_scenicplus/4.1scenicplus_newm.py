##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inferring eGRNs using SCENIC+: all NC cells
## Zhiyuan Hu
## 24 jan 2023
## last modified 15 dec 2023
##
## Env: scenicplus
## See README
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus")
# Refer to https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html
# /ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/4.1scenicplus_newm.py
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
work_dir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output'
tmpDir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/tmp'
old_work_dir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall/output'

# Read RNA data
rna = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall/scRNA/adata_ncall.h5ad'
adata = sc.read_h5ad(rna)
# adata.var_names = adata.var['_index']

# import pickle
# with open(work_dir + '../scRNA/metadata_bc.pkl', 'wb') as f:
#   pickle.dump(adata.obs_names, f)

# Read Cistopic obj and motif enrichment results
cistopic_obj = dill.load(open(os.path.join(work_dir, 'cisTopicObject.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adjust cell names
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

adata.obs['cell_type'] = adata.obs['cell_type'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
adata.obs['cell_type'].value_counts()
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace(' ', '_') for x in range(len(adata.obs['cell_type']))]
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace(':', '') for x in range(len(adata.obs['cell_type']))]
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace('.', '') for x in range(len(adata.obs['cell_type']))]
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace('+', '') for x in range(len(adata.obs['cell_type']))]
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace('-', '_') for x in range(len(adata.obs['cell_type']))]
adata.obs['cell_type'] = [adata.obs['cell_type'][x].replace('/', '_') for x in range(len(adata.obs['cell_type']))]

adata.obs['cell_type'].value_counts()


# Add barcode column
adata.obs['barcode'] = [x.split('_')[1] for x in adata.obs.index.tolist()]
adata.obs['sample_id'] = adata.obs['sample']
adata.obs['sample_id'] = adata.obs['sample_id'].replace('22cit', 's1')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('16dp', 's2')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('16cit', 's3')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('10dp', 's4')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('10cit', 's5')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('22dp', 's6')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('22cit2', 's7')
adata.obs['sample_id'] = adata.obs['sample_id'].replace('4mix', 's8')
adata.obs['sample_id'].value_counts()

# convert scRNA-seq barcodes to scATAC-seq ones
adata.obs_names = adata.obs['barcode'] + '___' + adata.obs['sample_id'].astype("object")

cistopic_obj.cell_names[0:5]
# adata2 = adata.X.to_adata()
# adata2.var_names = adata2.var['_index']
# 
# adata.layers['raw_counts']
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Create SCENIC+ object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Needs lots MEM(~80GB)  but the rstudio can still handle
# creating a scenicplus object containing all the analysis we have done up to this point
# load:
# 
# the AnnData object containing the scRNA-seq side of the analysis.
# the cisTopic object containing the scATAC-seq side of the analysis.
# the motif enrichment dictionary containing the motif enrichment results.

from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr
)
# 2023-01-24 06:43:55,986 cisTopic     INFO     Imputing drop-outs
# 2023-01-24 06:44:15,705 cisTopic     INFO     Scaling
# 2023-01-24 06:44:35,137 cisTopic     INFO     Keep non zero rows
# 2023-01-24 06:45:00,629 cisTopic     INFO     Imputed accessibility sparsity: 0.4743868317252342
# 2023-01-24 06:45:00,631 cisTopic     INFO     Create CistopicImputedFeatures object
# 2023-01-24 06:45:00,632 cisTopic     INFO     Done!
# del(adata)
del(cistopic_obj)

scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj

# #only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
# scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
# scplus_obj.dr_cell['GEX_rep'] = scplus_obj.dr_cell['GEX_X_tsne'].iloc[:, 0:2]

# filter low accessible regions and low expressed genes; the smallest population is 1.6% (C23)
filter_genes(scplus_obj, min_pct = 0.5)
# 2023-01-24 06:45:07,141 Preprocessing INFO     Going from 27599 genes to 17654 genes.
filter_regions(scplus_obj, min_pct = 0.5)
# 2023-01-24 06:51:40,665 Preprocessing INFO     Going from 374712 regions to 314399 regions.

print(scplus_obj)
# SCENIC+ object with n_cells x n_genes = 16550 x 17654 and n_cells x n_regions = 16550 x 314399

if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
    os.makedirs(os.path.join(work_dir, 'scenicplus'))

# Save
import pickle
with open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb') as f:
  pickle.dump(scplus_obj, f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.Generate cistromes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Merge cistromes (all)
from scenicplus.cistromes import *
import time
start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print("Time used to merge cistromes")
print(time/60) #11.719307907422383 minutes

# scplus_obj.uns['Cistromes']['Unfiltered'].keys()

# Save
import pickle
with open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb') as f:
  pickle.dump(scplus_obj, f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Infer enhancer to gene relationships
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

biomart_host = "http://dec2021.archive.ensembl.org/"

# download a list of known human TFs
# !wget -O pbmc_tutorial/data/utoronto_human_tfs_v_1.01.txt  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
# # download a the program bedToBigBed this will be used to generate files which can be uploaded to the UCSC genome browser
# !wget -O /ceph/home/z/zhu/t1data/multiome/R/GRN_SCENICplus/test/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
# !chmod +x /ceph/home/z/zhu/t1data/multiome/R/GRN_SCENICplus/test/bedToBigBed
 
############################################
# Get chromosome sizes (for danRer11 here) #
############################################
import pyranges as pr
import requests
# target_url='http://hgdownload.cse.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.chrom.sizes     ' # modified to zebrafish
target_url="/ceph/home/z/zhu/t1data/ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star/chrNameLength.txt"
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
# chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('_random', '') for x in range(len(chromsizes['Chromosome']))]
# chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
# chromsizes['Chromosome'] = [chromsizes['Chromosome'][x]+'.1' if 'chr' not in chromsizes['Chromosome'][x] else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)
chromsizes

#######################
# Get Gene annotations #
#######################
import pybiomart as pbm
# For zebrafish
dataset = pbm.Dataset(name='drerio_gene_ensembl',  host='http://dec2021.archive.ensembl.org/')
import pandas as pd
annot = pd.read_csv("/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105_scplus4.1_geneNameMatched.csv")
annot.Strand[annot.Strand == 1] = '+'
annot.Strand[annot.Strand == -1] = '-'
annot = pr.PyRanges(annot)
annot

# A. Get search space ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from scenicplus.enhancer_to_gene import calculate_regions_to_genes_relationships, GBM_KWARGS
sys.path.append('/ceph/home/z/zhu/t1data/multiome/R/GRN_SCENICplus/code')
import scenicplus_functions
scenicplus_functions.get_search_space_my(scplus_obj,
                 pr_annot=annot,
                 pr_chromsizes=chromsizes,
                 upstream = [1000, 150000],
                 downstream = [1000, 150000])
                 
scplus_obj.uns['search_space']

## Append foxd3 Cherry and foxd3 citrine
tmp=scplus_obj.uns['search_space'][scplus_obj.uns['search_space'].Gene == 'foxd3']
tmp.Gene='foxd3-mCherry'
scplus_obj.uns['search_space'] = scplus_obj.uns['search_space'].append(tmp)
tmp.Gene='foxd3-citrine'
scplus_obj.uns['search_space'] = scplus_obj.uns['search_space'].append(tmp)


if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
    os.makedirs(os.path.join(work_dir, 'scenicplus'))

# Save
import pickle
with open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb') as f:
  pickle.dump(scplus_obj, f)


# B. Enhancer-to-gene models ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Infer TF to gene relationships (Done by scripts)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Infer TF to gene relationships (Done by scripts)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

