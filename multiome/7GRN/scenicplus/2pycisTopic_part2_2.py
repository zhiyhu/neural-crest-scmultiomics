##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pycisTopic process for SCENIC+: NC all cells
## Zhiyuan Hu
## 16 Nov 2022
## last modified 23 oct 2023
## RUN by scripts/run_cistopic2.sh
## Env: scenicplus
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reticulate::use_condaenv("scenicplus")
# /ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/2pycisTopic_part2_2.py

# reference: https://pycistopic.readthedocs.io/en/latest/Cortex_pycisTopic.html#9.-Differentially-accessible-regions-(DARs)
# reference: https://scenicplus.readthedocs.io/en/latest/pbmc_multiome_tutorial.html#scATAC-seq-preprocessing-using-pycisTopic 

import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
pd.__version__
import pycisTopic
import matplotlib.pyplot as plt
pycisTopic.__version__
# Project directory
projDir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/'

# Output directory
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

# Temp dir
## this needs to be short to avoid "OSError: AF_UNIX path length cannot exceed 107 bytes"
tmpDir = '/ceph/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/tmp/'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. Topic modeling (Topic modeling can be computationaly intense!)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~ Run Model ~~~~~~~~~~~~~~~#

# Run by shell script (scripts/run_runModel.sh)

#~~~~~~~~~~ Model selection ~~~~~~~~~~#

# # Load cisTopic object
# import pickle
# infile = open(outDir + 'cisTopicObject.pkl', 'rb')
# cistopic_obj = pickle.load(infile)
# infile.close()
# 
# # Load models
# import pickle
# models = []
# for n_topics in [5, 10 ,15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
#     # infile = open(outDir+ 'models/mallet_nt' + str(n_topics) + '.pkl', 'rb')
#     infile = open(tmpDir + 'runModels/nt' + str(n_topics) + '/Topic' + str(n_topics) + '.pkl', 'rb')
#     models.append(pickle.load(infile)) 
#     infile.close()
# 
# if not os.path.exists(os.path.join(outDir,"models")):
#     os.makedirs(os.path.join(outDir,"models"))
# 
# from pycisTopic.lda_models import *
# model=evaluate_models(models,
#                      select_model=90,
#                      return_model=True,
#                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
#                      plot_metrics=False,
#                      save= outDir + 'models/model_selection.pdf')
# 
# # Add model to cisTopicObject
# cistopic_obj.add_LDA_model(model)
# 
# # Save
# with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
#   pickle.dump(cistopic_obj, f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 6. Visualization
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load cisTopic object
import pickle
infile = open(outDir + 'cisTopicObject.pkl', 'rb')
cistopic_obj = pickle.load(infile)
infile.close()

from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
run_tsne(cistopic_obj,target  = 'cell', scale=True)

if not os.path.exists(os.path.join(outDir, 'visualization')):
  os.mkdir(outDir+'/visualization')
  
plot_metadata(cistopic_obj,
                 reduction_name='UMAP',
                 variables=['cell_type'], # Labels from RNA and new clusters
                 target='cell', num_columns=1,
                 text_size=10,
                 dot_size=5,
                 figsize=(10,10),
                 save= outDir + 'visualization/dimensionality_reduction_umap_label.pdf')

plot_metadata(cistopic_obj,
                 reduction_name='tSNE',
                 variables=['cell_type','stage','genotype_new'], # Labels from RNA and new clusters
                 target='cell', num_columns=3,
                 text_size=10,
                 dot_size=5,
                 figsize=(30,10),
                 save= outDir + 'visualization/dimensionality_reduction_tsne_label.pdf')

# plot the topic-contributions
plot_topic(cistopic_obj,
            reduction_name = 'tSNE',
            target = 'cell',
            num_columns=5,
            save= outDir + 'visualization/dimensionality_reduction_tSNE_topic_contr.pdf')

# Heatmap
from pycisTopic.clust_vis import *
cell_topic_heatmap(cistopic_obj,
                     variables = ['cell_type'],
                     scale = False,
                     legend_loc_x = 1.05,
                     legend_loc_y = -1.2,
                     legend_dist_y = -1,
                     figsize=(10,20),
                     save = outDir + 'visualization/heatmap_topic_contr.pdf')

# Save
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 7. Inferring candidate enhancer regions (~2h)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7.1 binarization of region-topic probabilites.
if not os.path.exists(os.path.join(outDir, 'topic_binarization')):
    os.makedirs(os.path.join(outDir, 'topic_binarization'))

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu', plot=True, num_columns=5, save= outDir + 'topic_binarization/otsu.pdf')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

## To add: compute the topic quality control metrics????

# 7.2 calculation DARs per cell type.
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
len(variable_regions)
# 49812
## Slow step 
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, 
   variable='cell_type', 
   var_features=variable_regions, split_pattern = '-')

# Save
if not os.path.exists(os.path.join(outDir, 'candidate_enhancers')):
    os.makedirs(os.path.join(outDir, 'candidate_enhancers'))

import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(outDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(outDir, 'candidate_enhancers/markers_dict.pkl'), 'wb'))



