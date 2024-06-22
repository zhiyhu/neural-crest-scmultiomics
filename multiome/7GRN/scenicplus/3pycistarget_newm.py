##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pycisTarget process for SCENIC+
## Zhiyuan Hu
## 18 Nov 2022
## last modified 24 oct 2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reticulate::use_condaenv("scenicplus")
# multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/3pycistarget_newm.py
import pycistarget
import pandas as pd 
import os
pycistarget.__version__


# Project directory
projDir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/'

# Output directory
work_dir=projDir + 'output/'
proj='ncall_2023oct_ccb'

# tmp directort
tmpDir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/tmp/'
tmp_dir='multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/tmp/'
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1. Cistarget databases
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Local computer: preprocessing Danio-Code data
# prepare_cbuster_motifs.py
# 3create_cisTarget_databases.sh

# Define rankings, score and motif annotation database
db_fpath = work_dir + "ctx_db/"  # work_dir=projDir + 'output/'
motif_annot_fpath = "multiome/analysis_newref/GRN_scenicplus/data/motif_2023oct/"

rankings_db = db_fpath + 'ncall_newmotif' + '.regions_vs_motifs.rankings.feather'
scores_db =  db_fpath + 'ncall_newmotif' + '.regions_vs_motifs.scores.feather'
motif_annotation = os.path.join(motif_annot_fpath, 'dr11_motif2tf_curated_filtered2023oct.tbl')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2. Motif enrichment analysis 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pd.set_option('display.max_columns', 0)
motif2tf=pd.read_csv("multiome/analysis_newref/GRN_scenicplus/data/motif/dr11_motif2tf_curated_filtered20230122.tbl", sep = "\t")
tmp = pd.DataFrame([['glam2out_motif1','glam2out_motif1','glam2out_motif1', 'GLAM2', '1.0', 'foxd3','0',	'None',	'None',	'1',	'None',	'None',	'gene is directly annotated']], columns=motif2tf.columns)
motif2tf = motif2tf.append(tmp)
tmp = pd.DataFrame([['ets2','ets2','ets2', 'custom', '1.0', 'ets1','0',	'None',	'None',	'1',	'None',	'None',	'gene is directly annotated']], columns=motif2tf.columns)
motif2tf = motif2tf.append(tmp)
motif2tf 
motif2tf[898:899]
motif2tf.to_csv(motif_annotation, sep = "\t")

################################
# Load region binarized topics #
################################
import pickle
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'candidate_enhancers/markers_dict.pkl'), 'rb'))

#############################################
# Convert to dictionary of pyranges objects #
#############################################
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index #[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index #[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index #[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

#############################
#### Get TSS annotations ####
#############################
import pybiomart as pbm
# Annotation of ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
custom_annot = pd.read_csv("multiome/analysis_newref/GRN_scenicplus/data/gene_annot/pybiomart_drerio_gene_ensembl105_geneNameMatched.csv")
custom_annot

################################################################
#  run pycistarget using the run_pycistarget wrapper function  #
################################################################
if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))
from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    save_path = os.path.join(work_dir, 'motifs'),
    custom_annot = custom_annot,
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = '1.0'
    )

