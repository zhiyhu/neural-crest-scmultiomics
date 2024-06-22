##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pycisTopic process for SCENIC+: NC all cells
## Zhiyuan Hu
## 16 Nov 2022
## last modified 13 Dec 2023
##
## Env: scenicplus
## content: pseudobulk > consensus peaks > QC > cisTopic object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reticulate::use_condaenv("scenicplus_new")

# reference: https://pycistopic.readthedocs.io/en/latest/Cortex_pycisTopic.html#9.-Differentially-accessible-regions-(DARs)
# reference: https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#scATAC-seq-preprocessing-using-pycisTopic 

# multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/2pycisTopic_part1_1.py

import warnings
warnings.simplefilter(action='ignore')
import pandas as pd

import pycisTopic
import matplotlib.pyplot as plt
pycisTopic.__version__
# '1.0.2.dev21+g219225d'
import os 
# Project directories
projDir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/'
rna = os.path.join('multiome/analysis_newref/GRN_scenicplus/ncall/',"scRNA/adata_ncall.h5ad") # rna data
tmpDir = 'tmp/'## tmp dir needs to be short to avoid "OSError: AF_UNIX path length cannot exceed 107 bytes"

## Path to fragments files of samples
fragments_path = 'multiome/analysis_newref/cellranger_arc/output/'
fragments_dict = {'s1': fragments_path+'scmo_s1/outs/atac_fragments.tsv.gz',
                 's2': fragments_path+'scmo_s2/outs/atac_fragments.tsv.gz',
                 's3': fragments_path+'scmo_s3/outs/atac_fragments.tsv.gz',
                 's4': fragments_path+'scmo_s4/outs/atac_fragments.tsv.gz',
                 's5': fragments_path+'scmo_s5/outs/atac_fragments.tsv.gz',
                 's6': fragments_path+'scmo_s6/outs/atac_fragments.tsv.gz',
                 's7': fragments_path+'scmo_s7/outs/atac_fragments.tsv.gz',
                 's8': fragments_path+'scmo_s8/outs/atac_fragments.tsv.gz'}

# Output directory
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1. Getting pseudobulk profiles from cell annotations (~2h)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get metadata from a h5ad file
import scanpy as sc
adata = sc.read_h5ad(rna)
cell_data = adata.obs
del(adata)

cell_data['cell_type'] = cell_data['cell_type'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
cell_data['cell_type'].value_counts()

# tidy up cell type just in case
cell_data['cell_type'] = [cell_data['cell_type'][x].replace(' ', '_') for x in range(len(cell_data['cell_type']))]
cell_data['cell_type'] = [cell_data['cell_type'][x].replace(':', '') for x in range(len(cell_data['cell_type']))]
cell_data['cell_type'] = [cell_data['cell_type'][x].replace('.', '') for x in range(len(cell_data['cell_type']))]
cell_data['cell_type'] = [cell_data['cell_type'][x].replace('+', '') for x in range(len(cell_data['cell_type']))]
cell_data['cell_type'] = [cell_data['cell_type'][x].replace('-', '_') for x in range(len(cell_data['cell_type']))]
cell_data['cell_type'] = [cell_data['cell_type'][x].replace('/', '_') for x in range(len(cell_data['cell_type']))]

cell_data['cell_type'].value_counts()

# Add barcode column
cell_data['barcode'] = [x.split('_')[1] for x in cell_data.index.tolist()]
cell_data['sample_id'] = cell_data['sample']
cell_data['sample_id'] = cell_data['sample_id'].replace('scmo_s', 's')
cell_data['sample_id'] = cell_data['sample_id'].replace('22cit', 's1')
cell_data['sample_id'] = cell_data['sample_id'].replace('16dp', 's2')
cell_data['sample_id'] = cell_data['sample_id'].replace('16cit', 's3')
cell_data['sample_id'] = cell_data['sample_id'].replace('10dp', 's4')
cell_data['sample_id'] = cell_data['sample_id'].replace('10cit', 's5')
cell_data['sample_id'] = cell_data['sample_id'].replace('22dp', 's6')
cell_data['sample_id'] = cell_data['sample_id'].replace('22cit2', 's7')
cell_data['sample_id'] = cell_data['sample_id'].replace('4mix', 's8')
cell_data['sample_id'].value_counts()

# Get chromosome sizes (for danRer11 here)
import pyranges as pr
import requests
target_url="ref/ensembl105/cellranger_arc/GRCz11_ensembl105_foxd3_cellranger_arc/star/chrNameLength.txt"
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]

chromsizes=pr.PyRanges(chromsizes)
chromsizes
# +---------------+-----------+-----------+
# | Chromosome    | Start     | End       |
# | (category)    | (int64)   | (int64)   |
# |---------------+-----------+-----------|
# | 1             | 0         | 59578282  |
# | 2             | 0         | 59640629  |
# | 3             | 0         | 62628489  |
# | 4             | 0         | 78093715  |
# | ...           | ...       | ...       |
# | KZ116067.1    | 0         | 159186    |
# | MT            | 0         | 16596     |
# | foxd3-citrine | 0         | 1320      |
# | foxd3-mCherry | 0         | 1504      |
# +---------------+-----------+-----------+
# Unstranded PyRanges object has 995 rows and 3 columns from 995 chromosomes.
# For printing, the PyRanges was sorted on Chromosome.

# produce bigwig file
if not os.path.exists(outDir + 'consensus_peak_calling'):
    os.mkdir(outDir + 'consensus_peak_calling')
if not os.path.exists(outDir + 'consensus_peak_calling/pseudobulk_bed_files/'):
    os.mkdir(outDir + 'consensus_peak_calling/pseudobulk_bed_files/')
if not os.path.exists(outDir + 'consensus_peak_calling/pseudobulk_bw_files/'):
    os.mkdir(outDir + 'consensus_peak_calling/pseudobulk_bw_files/')

# this step takes ~ 30min
from pycisTopic.pseudobulk_peak_calling import *
bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                 variable = 'cell_type', # cell type annotation
                 sample_id_col = 'sample_id',
                 chromsizes = chromsizes,
                 bed_path = outDir + 'consensus_peak_calling/pseudobulk_bed_files/',
                 bigwig_path = outDir + 'consensus_peak_calling/pseudobulk_bw_files/',
                 path_to_fragments = fragments_dict,
                 n_cpu = 1,
                 normalize_bigwig = True,
                 remove_duplicates = True,
                 _temp_dir = tmpDir + 'ray_spill',
                 split_pattern = '___')
# Save
import pickle
with open(outDir + 'consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl', 'wb') as f:
  pickle.dump(bed_paths, f)

with open(outDir + 'consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl', 'wb') as f:
  pickle.dump(bw_paths, f)
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2. Inferring consensus peaks (~10min)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
from pycisTopic.pseudobulk_peak_calling import *
macs_path=".conda/envs/macs2_new/bin/macs2"
macs_outdir = outDir + 'consensus_peak_calling/MACS/'
if not os.path.exists(macs_outdir):
  os.mkdir(macs_outdir)
# Run peak calling ~10min
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_paths,
                                 macs_outdir,
                                 genome_size=1.4e9, ## modified this to zebrafish genome size
                                 n_cpu=1,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05,
                                 _temp_dir = tmpDir + 'ray_spill')
                                 
# Save
import pickle
with open(outDir + 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl', 'wb') as f:
  pickle.dump(narrow_peaks_dict, f)

from pycisTopic.iterative_peak_calling import *
# Other param
peak_half_width=250
path_to_blacklist="ref/blacklist/danRer11_blacklist_USCSliftover_20221223_noCHR.bed"

# Get consensus peaks
consensus_peaks=get_consensus_peaks(narrow_peaks_dict, peak_half_width, 
                                    chromsizes=chromsizes, 
                                    path_to_blacklist=path_to_blacklist)
consensus_peaks
# Write to bed
consensus_peaks.to_bed(path= outDir + 'consensus_peak_calling/consensus_regions.bed', keep=True, compression='infer', chain=False)



