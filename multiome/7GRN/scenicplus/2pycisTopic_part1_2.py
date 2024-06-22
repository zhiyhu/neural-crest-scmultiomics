##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pycisTopic process for SCENIC+: NC all cells
## Zhiyuan Hu
## 16 Nov 2022
## last modified 13 Dec 2023
##
## Env: scenicplus
## content: pseudobulk > consensus peaks > QC > cisTopic object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reticulate::use_condaenv("scenicplus")
# multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/2pycisTopic_part1_2.py

# reference: https://pycistopic.readthedocs.io/en/latest/Cortex_pycisTopic.html#9.-Differentially-accessible-regions-(DARs)
# reference: https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#scATAC-seq-preprocessing-using-pycisTopic 

import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import pycisTopic
import matplotlib.pyplot as plt
pycisTopic.__version__
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. QC (~3h)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sample level and cell-barcode level quality control statistics:
## Log number of unique fragments per cell barcode.
## FRIP per cell barcode.
## TSS enrichment per cell barcode.
## Duplication rate per cell barcode.

# Get TSS annotations
import pybiomart as pbm
# available genenames

import pandas as pd
genenames =  pd.read_csv("multiome/analysis_newref/preprocessing/rds/seu_featuredata.tsv", sep = '\t')

annot = pd.read_csv("multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105_geneNameMatched.csv")
# If you want to run all (or several of) the metrics, you can use the compute_qc_stats() function. 
# As input you need to provide a dictionary containing the fragments files per sample and another 
# dictionary the corresponding regions to use to estimate the FRIP.
from pycisTopic.qc import *
## Set regions. We will use the consensus peaks we have just called, but we could also use the bulk peaks per sample instead for this step
path_to_regions= {'s1': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's2': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's3': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's4': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's5': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's6': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's7': outDir + 'consensus_peak_calling/consensus_regions.bed',
                  's8': outDir + 'consensus_peak_calling/consensus_regions.bed',}
                  
# ~1h30min
metadata_bc, profile_data_dict = compute_qc_stats(fragments_dict = fragments_dict,
                tss_annotation = annot,
                stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
                label_list = None,
                path_to_regions = path_to_regions,
                n_cpu = 1,
                valid_bc = None,
                n_frag = 100,
                n_bc = None,
                tss_flank_window = 1000,
                tss_window = 50,
                tss_minimum_signal_window = 100,
                tss_rolling_window = 10,
                remove_duplicates = True,
                _temp_dir = tmpDir + 'ray_spill')

if not os.path.exists(os.path.join(outDir, 'quality_control')):
  os.makedirs(outDir+'quality_control')

import pickle
with open(outDir + 'quality_control/metadata_bc.pkl', 'wb') as f:
  pickle.dump(metadata_bc, f)

import pickle
with open(outDir + 'quality_control/profile_data_dict.pkl', 'wb') as f:
  pickle.dump(profile_data_dict, f)

#~~~~~~~~~~~~~~~~~~~~~~~~~ 3a. Sample-level statistics

# Load sample metrics
import pickle
infile = open(outDir + 'quality_control/profile_data_dict.pkl', 'rb')
profile_data_dict = pickle.load(infile)
infile.close()

from pycisTopic.qc import *
plot_sample_metrics(profile_data_dict,
           insert_size_distribution_xlim=[0,600],
           ncol=5,
           plot=True,
           save= outDir + 'quality_control/sample_metrics.pdf',
           duplicate_rate_as_hexbin = True)

#~~~~~~~~~~~~~~~~~~~~~~~~~ 3b. Barcode level statistics

# Load barcode metrics
import pickle
infile = open(outDir + 'quality_control/metadata_bc.pkl', 'rb')
metadata_bc = pickle.load(infile)
infile.close()

                         #[min,  #max]
QC_filters = {
    'Log_unique_nr_frag': [3 , None],
    'FRIP':               [0.35, None],
    'TSS_enrichment':     [4   , None],
    'Dupl_rate':          [None, None]

}

# Return figure to plot together with other metrics, and cells passing filters. Figure will be saved as pdf.
from pycisTopic.qc import *
FRIP_NR_FRAG_fig = {}
FRIP_NR_FRAG_filter = {}
TSS_NR_FRAG_fig = {}
TSS_NR_FRAG_filter = {}
DR_NR_FRAG_fig = {}
for sample in metadata_bc.keys():
    FRIP_NR_FRAG_fig[sample], FRIP_NR_FRAG_filter[sample]=plot_barcode_metrics(metadata_bc[sample],
                                           var_x='Log_unique_nr_frag',
                                           var_y='FRIP',
                                           min_x=QC_filters['Log_unique_nr_frag'][0],
                                           max_x=QC_filters['Log_unique_nr_frag'][1],
                                           min_y=QC_filters['FRIP'][0],
                                           max_y=QC_filters['FRIP'][1],
                                           return_cells=True,
                                           return_fig=True,
                                           plot=False,
                                           save= outDir + 'quality_control/barcode_metrics_FRIP-VS-NRFRAG_'+sample+'.pdf')
    # Return figure to plot together with other metrics, and cells passing filters
    TSS_NR_FRAG_fig[sample], TSS_NR_FRAG_filter[sample]=plot_barcode_metrics(metadata_bc[sample],
                                          var_x='Log_unique_nr_frag',
                                          var_y='TSS_enrichment',
                                          min_x=QC_filters['Log_unique_nr_frag'][0],
                                          max_x=QC_filters['Log_unique_nr_frag'][1],
                                          min_y=QC_filters['TSS_enrichment'][0],
                                          max_y=QC_filters['TSS_enrichment'][1],
                                          return_cells=True,
                                          return_fig=True,
                                          plot=False,
                                          save= outDir + 'quality_control/barcode_metrics_TSS-VS-NRFRAG_'+sample+'.pdf')
    # Return figure to plot together with other metrics, but not returning cells (no filter applied for the duplication rate  per barcode)
    DR_NR_FRAG_fig[sample]=plot_barcode_metrics(metadata_bc[sample],
                                          var_x='Log_unique_nr_frag',
                                          var_y='Dupl_rate',
                                          min_x=QC_filters['Log_unique_nr_frag'][0],
                                          max_x=QC_filters['Log_unique_nr_frag'][1],
                                          min_y=QC_filters['Dupl_rate'][0],
                                          max_y=QC_filters['Dupl_rate'][1],
                                          return_cells=False,
                                          return_fig=True,
                                          plot=False,
                                          plot_as_hexbin = True)


# Plot barcode stats in one figure
fig=plt.figure(figsize=(40, 80))
i=1
for sample in FRIP_NR_FRAG_fig.keys():
    plt.subplot(8, 3, i)
    plt.gca().set_title(sample, fontsize=20)
    i += 1
    img = fig2img(FRIP_NR_FRAG_fig[sample]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(8, 3, i)
    plt.gca().set_title(sample, fontsize=20)
    i += 1
    img = fig2img(TSS_NR_FRAG_fig[sample])
    plt.imshow(img)
    plt.axis('off')
    plt.subplot(8, 3, i)
    plt.gca().set_title(sample, fontsize=20)
    i += 1
    img = fig2img(DR_NR_FRAG_fig[sample])
    plt.imshow(img)
    plt.axis('off')
plt.savefig(outDir + 'quality_control/combined_qc.pdf')

# select the cells passing filters
sel_cells_dict = {}
for sample in FRIP_NR_FRAG_filter.keys():
    sel_cells_dict[sample] = list((set(FRIP_NR_FRAG_filter[sample]) & set(TSS_NR_FRAG_filter[sample])))

# check how many match with the scRNA-seq part of the multiome
cell_data_small = cell_data[['sample_id', 'barcode']]
sel_cells_rna = {x: cell_data_small[cell_data_small['sample_id'] == x]['barcode'].tolist() for x in set(cell_data_small['sample_id'])}

for sample in sel_cells_dict.keys():
    sel_cells_dict[sample] = list((set(sel_cells_dict[sample]) & set(sel_cells_rna[sample])))



# keep cells that have high quality profiles in both the scATAC-seq and scRNA-seq data
import pickle
with open(outDir +'/quality_control/bc_passing_filters.pkl', 'wb') as f:
  pickle.dump(sel_cells_dict, f)

import pickle
infile = open(outDir +'/quality_control/bc_passing_filters.pkl', 'rb')
bc_passing_filters = pickle.load(infile)
infile.close()

for sample in sel_cells_dict.keys():
    print(f"{sample}: {len(bc_passing_filters[sample])} barcodes passed QC stats")

sum_bc=0
for sample in sel_cells_dict.keys():
    sum_bc=sum_bc+len(bc_passing_filters[sample])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 4. Creating a cisTopic object 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate a binary count matrix of ATAC-seq fragments over consensus peaks

# Path to regions
path_to_regions = outDir + 'consensus_peak_calling/consensus_regions.bed'
## Blacklist
path_to_blacklist="ref/blacklist/danRer11_blacklist_USCSliftover_20221223_noCHR.bed"
# Metrics
import pickle
infile = open(outDir + 'quality_control/metadata_bc.pkl', 'rb')
metadata_bc = pickle.load(infile)
infile.close()
# Valid barcodes
import pickle
infile = open(outDir +'/quality_control/bc_passing_filters.pkl', 'rb')
bc_passing_filters = pickle.load(infile)
infile.close()
# Create cisTopic object
from pycisTopic.cistopic_class import *
# takes ~1h
cistopic_obj_list=[create_cistopic_object_from_fragments(path_to_fragments=fragments_dict[key],
                                               path_to_regions=path_to_regions,
                                               path_to_blacklist=path_to_blacklist,
                                               metrics=metadata_bc[key],
                                               valid_bc=bc_passing_filters[key],
                                               n_cpu=1,
                                               project=key) for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)

# Cell metadata
cell_data.index = cell_data['barcode'] + '___' + cell_data['sample_id'].astype("object")
cistopic_obj.add_cell_data(cell_data)

print(cistopic_obj)

# Save
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)
