##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Explore GRN
## Zhiyuan Hu
## 28 Nov 2022
## last modified 20 dec 2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus_new")
# Refer to https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html
# /ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/5expore_grn.py

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
work_dir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/output/'
tmpDir = '/ceph/home/z/zhu/t1data/tmp'

import os
if not os.path.exists(os.path.join(work_dir, 'scenicplus/plots')):
    os.makedirs(os.path.join(work_dir, 'scenicplus/plots'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. Exploring SCENIC+ results
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pickle
infile = open(os.path.join(work_dir, 'scenicplus/scplus_obj_grn.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A. Generate eRegulon metadata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', 
    TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')

import pandas as pd
pd.set_option('display.max_columns', None)
scplus_obj.uns['eRegulon_metadata'][0:10]
scplus_obj.uns['eRegulon_metadata']["Gene_signature_name"].value_counts()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B. Assesing eGRN enrichment in cells
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format eRegulons
from scenicplus.eregulon_enrichment import *
get_eRegulons_as_signatures(scplus_obj, eRegulon_metadata_key='eRegulon_metadata', key_added='eRegulon_signatures')

## Score chromatin layer
# Region based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
region_ranking = make_rankings(scplus_obj, target='region')
# Score region regulons
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 1)
time = time.time()-start_time
print(time/60)
# 11.083547727266948

## Score transcriptome layer
# Gene based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
gene_ranking = make_rankings(scplus_obj, target='gene')
# Score gene regulons
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 1)
time = time.time()-start_time
print(time/60)
# 0.28899283011754356

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C. Assessing TF-eGRN relationships
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate pseudobulks
import time
start_time = time.time()
generate_pseudobulks(scplus_obj,
                         variable = 'ACC_cell_type',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Gene_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)
generate_pseudobulks(scplus_obj,
                         variable = 'ACC_cell_type',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Region_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)
time = time.time()-start_time
print(time/60)
# 0.6083502411842346

# Correlation between TF and eRegulons
import time
start_time = time.time()
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_cell_type',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Gene_based',
                        out_key = 'ACC_cell_type_eGRN_gene_based')
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_cell_type',
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Region_based',
                        out_key = 'ACC_cell_type_eGRN_region_based')
time = time.time()-start_time
print(time/60)
# 0.012485591570536296

scplus_obj.uns['Pseudobulk']['ACC_cell_type']['eRegulon_AUC']['Gene_based'][0:10]

# plot eRegulon enrichment versus TF expression for each pseudobulk
# Region based
# %matplotlib inline
import seaborn as sns
sns.set_style("white")
fig, ax = plt.subplots(1,1, figsize=(6,6))
prune_plot(scplus_obj,
           name = 'mitfa_+_+',
           pseudobulk_variable = 'ACC_cell_type',
           show_dot_plot = True,
           show_line_plot = False,
           # color_dict = color_dict,
           use_pseudobulk = True,
           auc_key = 'eRegulon_AUC',
           signature_key = 'Region_based',
           seed=555, ax = ax)
fig.tight_layout()
plt.show()
plt.savefig(work_dir + 'scenicplus/plots/eRegulon_enrichment_versus_TFexpr_regionbased_mitfa.pdf')

# Gene based
# %matplotlib inline
sns.set_style("white")
# colors = ["#E9842C","#F8766D", "#BC9D00", "#00C0B4", "#9CA700", "#6FB000", "#00B813", "#00BD61", "#00C08E", "#00BDD4",
#            "#00A7FF", "#7F96FF", "#E26EF7", "#FF62BF", "#D69100", "#BC81FF"]
# categories = sorted(set(scplus_obj.metadata_cell['ACC_cell_type']))
# color_dict = dict(zip(categories, colors[0:len(categories)]))

fig, ax = plt.subplots(1,1, figsize=(6,6))
prune_plot(scplus_obj,
           'mitfa_+_+',
           pseudobulk_variable = 'ACC_cell_type',
           show_dot_plot = True,
           show_line_plot = False,
           # color_dict = color_dict,
           use_pseudobulk = True,
           auc_key = 'eRegulon_AUC',
           signature_key = 'Gene_based',
           seed=555 , ax = ax)
fig.tight_layout()
plt.show()
plt.savefig(work_dir + 'scenicplus/plots/eRegulon_enrichment_versus_TFexpr_genebased_mitfa.pdf')


fig, ax = plt.subplots(1,1, figsize=(6,6))
prune_plot(scplus_obj,
           'sox6_+_+',
           pseudobulk_variable = 'ACC_cell_type',
           show_dot_plot = True,
           show_line_plot = False,
           # color_dict = color_dict,
           use_pseudobulk = True,
           auc_key = 'eRegulon_AUC',
           signature_key = 'Gene_based',
           seed=555 , ax = ax)
fig.tight_layout()
plt.show()
plt.savefig(work_dir + 'scenicplus/plots/eRegulon_enrichment_versus_TFexpr_genebased_sox6.pdf')

fig, ax = plt.subplots(1,1, figsize=(6,6))
prune_plot(scplus_obj,
           'sox6_+_+',
           pseudobulk_variable = 'ACC_cell_type',
           show_dot_plot = True,
           show_line_plot = False,
           # color_dict = color_dict,
           use_pseudobulk = True,
           auc_key = 'eRegulon_AUC',
           signature_key = 'Region_based',
           seed=555 , ax = ax)
fig.tight_layout()
plt.show()
plt.savefig(work_dir + 'scenicplus/plots/eRegulon_enrichment_versus_TFexpr_regionbased_sox6.pdf')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D. Identification of high quality regulons
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Correlation between region based regulons and gene based regulons
print("D. Identification of high quality regulons")
import pandas
df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)
len(correlations)
# 822

# plot
plt.figure(0)
plt.hist(correlations, 20, density=True, facecolor='g', alpha=0.75)
plt.xlabel("Correlation between Gene_based and Region_based eRegulon_AUCs")
plt.title("Distribution")
plt.show()
plt.savefig(work_dir + 'scenicplus/plots/eRegulon_AUC_correlation.pdf')

correlations.to_csv(work_dir + 'scenicplus/plots/eRegulon_AUC_correlation.csv')
correlations = correlations[abs(correlations) > 0.4] ## modify this to 0.5
# 559

# Keep only R2G +
keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x]
# Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended
# Keep regulons with more than 5 genes
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 5)]
keep_all = [x.split('_(')[0] for x in keep_gene]
keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

len(keep_gene)
# 186
len(scplus_obj.uns['selected_eRegulons']['Gene_based'])
# 186
len(scplus_obj.uns['selected_eRegulons']['Region_based'])
# 186

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E. Overlap between eRegulons
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# import matplotlib
# matplotlib.pyplot.margins(0, 0) 

print("E. Overlap between eRegulons")

# assess which eRegulons tend to be enriched in the same group of cells
from scenicplus.plotting.correlation_plot import *
correlation_heatmap(scplus_obj,
                    auc_key = 'eRegulon_AUC',
                    signature_keys = ['Gene_based'],
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
                    fcluster_threshold = 0.1, figsize = (40,40),
                    fontsize = 2, use_plotly = False, save = work_dir + 'scenicplus/plots/correlation_heatmap.pdf')

# check the overlap between eRegulons
#from scenicplus.plotting.correlation_plot import *
mat_jaccard=jaccard_heatmap(scplus_obj,
                    gene_or_region_based = 'Gene_based',
                    signature_key = 'eRegulon_signatures',
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 2, figsize = (40,40), return_data=True,
                    method='intersect', use_plotly = False, save = work_dir + 'scenicplus/plots/jaccard_heatmap_gene_based.pdf')
mat_jaccard[0].to_csv(work_dir + "scenicplus/plots/jaccard_matrix_gene_based.csv")

scplus_obj.uns['selected_eRegulons']['Region_based']

mat_jaccard=jaccard_heatmap(scplus_obj,
                    gene_or_region_based = 'Region_based',
                    signature_key = 'eRegulon_signatures',
                    selected_regulons = scplus_obj.uns['selected_eRegulons']['Region_based'],
                    fcluster_threshold = 0.1,
                    fontsize = 2, figsize = (40,40), return_data=True,
                    method='intersect', use_plotly = False, save = work_dir + 'scenicplus/plots/jaccard_heatmap_region_based.pdf')
mat_jaccard[0].to_csv(work_dir + "scenicplus/plots/jaccard_matrix_region_based.csv")

# binarize the eRegulons as in SCENIC. This information will be used afterwards for generating the loom file.
# binarize_AUC(scplus_obj,
#              auc_key='eRegulon_AUC',
#              out_key='eRegulon_AUC_thresholds',
#              signature_keys=['Gene_based', 'Region_based'],
#              n_cpu=2)

# import dill
# with open(work_dir + 'scenicplus/scplus_obj_grn.pkl', 'wb') as f:
#   dill.dump(scplus_obj, f)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# F. eGRN dimensionality reduction
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("F. eGRN dimensionality reduction")
# import pickle
# infile = open(os.path.join(work_dir, 'scenicplus/scplus_obj_grn.pkl'), 'rb')
# scplus_obj = pickle.load(infile)
# infile.close()


# combination of both gene and region based eRegulons results in a better clustering.
from scenicplus.dimensionality_reduction import *
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Gene_based', 'Region_based'], 
                   reduction_name='eRegulons_UMAP',
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
run_eRegulons_tsne(scplus_obj,
                   scale=True, signature_keys=['Gene_based', 'Region_based'], 
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
# using the layers independently as well
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Gene_based'],
                   reduction_name='eRegulons_UMAP_gb', 
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
run_eRegulons_tsne(scplus_obj,
                   scale=True, signature_keys=['Gene_based'],
                   reduction_name='eRegulons_tSNE_gb', 
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'])
run_eRegulons_umap(scplus_obj,
                   scale=True, signature_keys=['Region_based'],
                   reduction_name='eRegulons_UMAP_rb', 
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Region_based'])
run_eRegulons_tsne(scplus_obj,
                   scale=True, signature_keys=['Region_based'],
                   reduction_name='eRegulons_tSNE_rb', 
                   selected_regulons=scplus_obj.uns['selected_eRegulons']['Region_based'])

plot_metadata(scplus_obj,
                 reduction_name='eRegulons_UMAP_rb',
                 variables=['ACC_cell_type'],
                 num_columns=1,
                 text_size=5,
                 dot_size=2, save = work_dir + 'scenicplus/plots/eRegulons_UMAP_rb_cell_type.pdf')

plot_metadata(scplus_obj,
                 reduction_name='eRegulons_tSNE_gb',
                 variables=['ACC_cell_type'],
                 num_columns=1,
                 text_size=5,
                 dot_size=2, save = work_dir + 'scenicplus/plots/eRegulons_tSNE_gb_cell_type.pdf')

plot_metadata(scplus_obj,
                 reduction_name='eRegulons_tSNE_rb',
                 variables=['ACC_cell_type'],
                 num_columns=1,
                 text_size=5,
                 dot_size=2, save = work_dir + 'scenicplus/plots/eRegulons_tSNE_rb_cell_type.pdf')

plot_metadata(scplus_obj,
                 reduction_name='eRegulons_UMAP_gb',
                 variables=['ACC_cell_type'],
                 num_columns=1,
                 text_size=5,
                 dot_size=2, save = work_dir + 'scenicplus/plots/eRegulons_UMAP_gb_cell_type.pdf')


plot_metadata(scplus_obj,
                 reduction_name='eRegulons_UMAP_gb',
                 variables=['ACC_cell_type'],
                 num_columns=1,
                 dot_size=2)



## plot eRegulons
# For example, for some OL TFs
plot_eRegulon(scplus_obj,
              reduction_name='eRegulons_tSNE_rb',
              selected_regulons=['sox6_+_+','sox9b_+_+', 'zeb1b_+_+'],
              normalize_tf_expression=True,
              dot_size = 5, save = work_dir + 'scenicplus/plots/eRegulons_tSNE_rb_sox6_sox9b_zeb1b.pdf')

plot_eRegulon(scplus_obj,
              reduction_name='eRegulons_tSNE_rb',
              selected_regulons=['sox10_+_+','tfec_+_+', 'tfeb_+_+', 'mitfa_+_+', 'bhlhe40_+_+'],
              normalize_tf_expression=True,
              dot_size = 5, save = work_dir + 'scenicplus/plots/eRegulons_tSNE_rb_pigment.pdf')


plot_eRegulon(scplus_obj,
              reduction_name='eRegulons_tSNE_rb',
              selected_regulons=['foxd3_+_+','foxd3_-_+', 'sox10_+_+', 'sox9b_+_+', 'sox4a_+_+', 'fli1a_+_+', 'twist1a_+_-'],
              normalize_tf_expression=True,
              dot_size = 2, save = work_dir + 'scenicplus/plots/eRegulons_tSNE_rb_foxd3.png')


import dill
with open(work_dir + 'scenicplus/scplus_obj_grn.pkl', 'wb') as f:
  dill.dump(scplus_obj, f)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# G. eRegulon specificity scores
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("G. eRegulon specificity scores")
from scenicplus.RSS import *
regulon_specificity_scores(scplus_obj,
                         'ACC_cell_type',
                         signature_keys=['Gene_based'],
                         selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'],
                         out_key_suffix='_gene_based',
                         scale=False)
plot_rss(scplus_obj, 'ACC_cell_type_gene_based', num_columns=4, top_n=10,
save = work_dir + 'scenicplus/plots/regulon_specificity_scores.pdf')

# write RSS to csv
scplus_obj.uns["RSS"]['ACC_cell_type_gene_based'].to_csv(work_dir + 'scenicplus/plots/regulon_specificity_scores_gene_based.csv')
scplus_obj.uns["RSS"]['ACC_cell_type_region_based'].to_csv(work_dir + 'scenicplus/plots/regulon_specificity_scores_region_based.csv')


from scenicplus.RSS import *
regulon_specificity_scores(scplus_obj,
                         'ACC_cell_type',
                         signature_keys=['Gene_based'],
                         selected_regulons=scplus_obj.uns['selected_eRegulons']['Gene_based'],
                         out_key_suffix='_gene_based',
                         scale=False)

# from scenicplus.plotting.dotplot import *
sys.path.append('/ceph/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/code/custom_functions/')
import dotplot_custom
dotplot_custom.heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.to_df('EXP'),
        color_matrix = scplus_obj.uns['RSS']['ACC_cell_type_gene_based'],
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'ACC_cell_type',
        subset_eRegulons = scplus_obj.uns['selected_eRegulons']['Gene_based'],
        figsize = (45, 7),
        orientation = 'horizontal',
        split_repressor_activator=True, 
        save = work_dir + 'scenicplus/plots/specificity_genebased.pdf')

from scenicplus.RSS import *
regulon_specificity_scores(scplus_obj,
                         'ACC_cell_type',
                         signature_keys=['Region_based'],
                         selected_regulons=scplus_obj.uns['selected_eRegulons']['Region_based'],
                         out_key_suffix='_region_based',
                         scale=False)

# from scenicplus.plotting.dotplot import *
dotplot_custom.heatmap_dotplot(
        scplus_obj = scplus_obj,
        size_matrix = scplus_obj.to_df('EXP'),
        color_matrix = scplus_obj.uns['RSS']['ACC_cell_type_region_based'],
        scale_size_matrix = True,
        scale_color_matrix = True,
        group_variable = 'ACC_cell_type',
        subset_eRegulons = scplus_obj.uns['selected_eRegulons']['Region_based'],
        figsize = (35, 7),
        orientation = 'horizontal',
        split_repressor_activator=False, 
        save = work_dir + 'scenicplus/plots/specificity_regionbased.pdf')

import dill
with open(work_dir + 'scenicplus/scplus_obj_grn.pkl', 'wb') as f:
  dill.dump(scplus_obj, f)

scplus_obj.to_df('EXP').to_csv( work_dir + 'scenicplus/plots/scplus_obj_exp.txt.gz', sep = '\t')
scplus_obj.uns['RSS']['ACC_cell_type_gene_based'].to_csv( work_dir + 'scenicplus/plots/scplus_obj_RSS_cell_type_gene_based.csv')
scplus_obj.uns['RSS']['ACC_cell_type_region_based'].to_csv( work_dir + 'scenicplus/plots/scplus_obj_RSS_cell_type_region_based.csv')

# open file in write mode
with open(work_dir + 'scenicplus/plots/scplus_obj_selected_eRegulons_gene_based.csv', 'w') as fp:
    for item in scplus_obj.uns['selected_eRegulons']['Gene_based']:
        # write each item on a new line
        fp.write("%s\n" % item)
    print('Done')

with open(work_dir + 'scenicplus/plots/scplus_obj_selected_eRegulons_region_based.csv', 'w') as fp:
    for item in scplus_obj.uns['selected_eRegulons']['Region_based']:
        # write each item on a new line
        fp.write("%s\n" % item)
    print('Done')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# H. Integrated multiome plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Generate interaction and annotation pyranges
# import matplotlib.pyplot as plt
# import os
# from scenicplus.utils import get_interaction_pr
# import pyranges as pr
# bigwig_dir = '/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/10x_multiome_brain/output/atac/pycistopic/consensus_peak_calling/seurat_pseudobulk_bw_files/'
# bw_dict = {x.replace('.bw', ''): os.path.join(bigwig_dir, x) for x in os.listdir(bigwig_dir) if '.bw' in x}
# pr_consensus_bed = pr.read_bed('/staging/leuven/stg_00002/lcb/cbravo/Multiomics_pipeline/analysis/10x_multiome_brain/output/atac/pycistopic/consensus_peak_calling/consensus_regions.bed')
# pr_interact = get_interaction_pr(scplus_obj, 'hsapiens', 'hg38', inplace = False, subset_for_eRegulons_regions = True, eRegulons_key = 'eRegulons_importance')
# gtf_file = "/lustre1/project/stg_00002/lcb/fderop/data/00000000_genomes/GRCh38_STAR_2.7.5_rna/genes.gtf"
# pr_gtf = pr.read_gtf(gtf_file)
# 
# # Plot
# from importlib import reload
# from scenicplus.plotting.coverageplot import *
# fig = coverage_plot(
#         SCENICPLUS_obj = scplus_obj,
#         bw_dict = bw_dict,
#         region = 'chr1:225800009-225860175',
#         figsize = (10,20),
#             pr_gtf = pr_gtf,
#         color_dict = None,
#         plot_order = None,
#         pr_interact = pr_interact,
#         genes_violin_plot = ['TMEM63A', 'SOX5'],
#         meta_data_key = 'ACC_Seurat_cell_type',
#         pr_consensus_bed = pr_consensus_bed,
#         fontsize_dict={'bigwig_label': 12, 'gene_label': 0, 'violinplots_xlabel': 10, 'title': 12, 'bigwig_tick_label': 0, 'violinplots_ylabel': 3},
#         height_ratios_dict = {'bigwig_violin': 1, 'genes': 0.5, 'arcs': 10, 'custom_ax': 5})
# plt.tight_layout()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. Adding DEGs and DARs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# from pyscenic.diff_features import *
# get_differential_features(scplus_obj, 'ACC_cell_type', use_hvg = True, contrast_type = ['DARs', 'DEGs'])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# J. Export to loom
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# L. Plotting networks
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# from scenicplus.networks import *
# import networkx as nx
# # subset_genes = [ 'fhod1', 'quo', 'plxnb2a', 'rps15a', 'arpc1a', 'rap1b', 'tead3a', 'bach2b', 'wnt5b', 'myo1ea', 'hmga2', 'rap1b', 'ctdspl2b', 'cdc42ep4a', 'glis2b', 'elavl3', 'tjp1a', 'wwc3', 'ERC1', 'flncb']
# 
# submetadata = scplus_obj.uns['eRegulon_metadata'][scplus_obj.uns['eRegulon_metadata']['Region_signature_name'] == 'sox10_+_+_(104r)']
# submetadata = scplus_obj.uns['eRegulon_metadata'][scplus_obj.uns['eRegulon_metadata']['Region_signature_name'] == 'mitfa_+_+_(87r)']
# 
# 
# subset_genes = set(subset_genes)
# # len(subset_genes)
# 
# subset_genes =  ['bhlhe40', 'fgfr4', 'ctcf', 'tfec','ptpn4a','erbb3b','celf2','cax1', 'sox10','mitfa', 'sptlc2b']
# 
# # celf2: nervous system
# # sptlc2b: sensory and autonomic neuropathy type 1C?
# # 
# 
# nx_tables = create_nx_tables(scplus_obj,
#                     eRegulon_metadata_key = 'eRegulon_metadata',
#                     subset_eRegulons = ['bhlhe40','tfec','sox10','mitfa'],
#                     subset_regions = None,
#                     subset_genes = subset_genes,
#                     add_differential_gene_expression = True,
#                     add_differential_region_accessibility = True,
#                     differential_variable = ['ACC_cell_type'])
# 
# from scenicplus.networks import *
# G_c, pos_c, edge_tables_c, node_tables_c = create_nx_graph(nx_tables,
#                    use_edge_tables = ['TF2R','R2G'],
#                    color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {'bhlhe40': 'Orange', 'tfec': 'Purple', 'sox10': "Blue", 'mitfa': "Pink"}},
#                                     'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
#                    transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
#                    width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
#                    color_node_by = {'TF': {'variable': 'TF', 'category_color' : {'bhlhe40': 'Orange', 'tfec': 'Purple', 'sox10': "Blue", 'mitfa': "Pink"}},
#                                     'Gene': {'variable': 'ACC_cell_type_Log2FC_0NC_M_mid', 'continuous_color' : 'bwr'},
#                                     'Region': {'variable': 'ACC_cell_type_Log2FC_0NC_M_mid', 'continuous_color' : 'viridis'}},
#                    transparency_node_by =  {'Region': {'variable' : 'ACC_cell_type_Log2FC_0NC_M_mid', 'min_alpha': 0.1},
#                                     'Gene': {'variable' : 'ACC_cell_type_Log2FC_0NC_M_mid', 'min_alpha': 0.1}},
#                    size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
#                                     'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
#                                     'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
#                    shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                     'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                     'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
#                    label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
#                                     'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
#                                     'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
#                    layout='concentrical_layout',
#                    scale_position_by=250)
# 
# # To draw with networkx
# plt.figure(figsize=(10,10))
# plot_networkx(G_c, pos_c)                 
# plt.savefig(work_dir+ "scenicplus/networkx/"+ "bhlhe40_tfec_sox10.pdf")                    
# 
# # Kamada Kawai algorithm
# G_kk, pos_kk, edge_tables_kk, node_tables_kk = create_nx_graph(nx_tables,
#                    use_edge_tables = ['TF2R','R2G'],
#                    color_edge_by = {'TF2R': {'variable' : 'TF', 'category_color' : {'bhlhe40': 'Orange', 'tfec': 'Purple', 'sox10': "Blue", 'mitfa': "Pink"}},
#                                     'R2G': {'variable' : 'R2G_rho', 'continuous_color' : 'viridis', 'v_min': -1, 'v_max': 1}},
#                    transparency_edge_by =  {'R2G': {'variable' : 'R2G_importance', 'min_alpha': 0.1, 'v_min': 0}},
#                    width_edge_by = {'R2G': {'variable' : 'R2G_importance', 'max_size' :  1.5, 'min_size' : 1}},
#                    color_node_by = {'TF': {'variable': 'TF', 'category_color' : {'bhlhe40': 'Orange', 'tfec': 'Purple', 'sox10': "Blue", 'mitfa': "Pink"}},
#                                     'Gene': {'variable': 'ACC_cell_type_Log2FC_0NC_M_mid', 'continuous_color' : 'bwr'},
#                                     'Region': {'variable': 'ACC_cell_type_Log2FC_0NC_M_mid', 'continuous_color' : 'viridis'}},
#                    transparency_node_by =  {'Region': {'variable' : 'ACC_cell_type_Log2FC_0NC_M_mid', 'min_alpha': 0.1},
#                                     'Gene': {'variable' : 'ACC_cell_type_Log2FC_0NC_M_mid', 'min_alpha': 0.1}},
#                    size_node_by = {'TF': {'variable': 'fixed_size', 'fixed_size': 30},
#                                     'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
#                                     'Region': {'variable': 'fixed_size', 'fixed_size': 10}},
#                    shape_node_by = {'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                     'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                     'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
#                    label_size_by = {'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
#                                     'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 10.0},
#                                     'Region': {'variable': 'fixed_label_size', 'fixed_label_size': 0.0}},
#                    layout='kamada_kawai_layout',
#                    scale_position_by = 500)
# 
# plt.figure(figsize=(10,10))
# plot_networkx(G_kk, pos_kk)
# 
# export_to_cytoscape(G_kk, pos_kk,
#                     work_dir+ "scenicplus/cytoscape/"+ "bhlhe40_tfec_sox10_kk.cyjs",
#                     pos_scaling_factor=200, size_scaling_factor=1)
# 
# 



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# M. Exporting to UCSC
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# N. Save the eRegulon data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
if not os.path.exists(work_dir + "scenicplus/cytoscape"):
    os.makedirs(work_dir + "scenicplus/cytoscape")
    
scplus_obj.uns['eRegulon_metadata'].to_csv(work_dir + "scenicplus/cytoscape/eRegulon_metadata_all.csv")
regulon_metadata = scplus_obj.uns['eRegulon_metadata']

regulon_metadata = regulon_metadata[regulon_metadata['Gene_signature_name'].isin(scplus_obj.uns['selected_eRegulons']['Gene_based'])]
regulon_metadata.to_csv(work_dir + "scenicplus/cytoscape/eRegulon_metadata_filtered.csv")

import dill
with open(work_dir + 'scenicplus/scplus_obj_grn.pkl', 'wb') as f:
  dill.dump(scplus_obj, f)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save AUCell
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not os.path.exists(work_dir + "scenicplus/auc"):
    os.makedirs(work_dir + "scenicplus/auc")
    
scplus_obj.uns['eRegulon_AUC']['Region_based'].to_csv(work_dir + 'scenicplus/auc/eRegulon_AUC_region_based.tsv.gz', sep = "\t")
scplus_obj.uns['eRegulon_AUC']['Gene_based'].to_csv(work_dir + 'scenicplus/auc/eRegulon_AUC_gene_based.tsv.gz', sep = "\t")

