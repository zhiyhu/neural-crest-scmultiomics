##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pycisTopic process for SCENIC+: NC all cells
## Zhiyuan Hu
## 16 Nov 2022
## last modified 14 Dec 2023
##
## Env: scenicplus
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reticulate::use_condaenv("scenicplus")

# multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/2pycisTopic_part2_1.py

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
projDir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/'

# Output directory
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

# Temp dir
## this needs to be short to avoid "OSError: AF_UNIX path length cannot exceed 107 bytes"
tmpDir = 'multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/tmp/'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. Topic modeling (Topic modeling can be computationaly intense!)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~ Run Model ~~~~~~~~~~~~~~~#

# Run by shell script (scripts/run_runModel.sh)

#~~~~~~~~~~ Model selection ~~~~~~~~~~#

# Load cisTopic object
import pickle
infile = open(outDir + 'cisTopicObject.pkl', 'rb')
cistopic_obj = pickle.load(infile)
infile.close()

# Load models
import pickle
models = []
for n_topics in [5, 10 ,15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]:
    # infile = open(outDir+ 'models/mallet_nt' + str(n_topics) + '.pkl', 'rb')
    infile = open(tmpDir + 'runModels/nt' + str(n_topics) + '/Topic' + str(n_topics) + '.pkl', 'rb')
    models.append(pickle.load(infile)) 
    infile.close()

if not os.path.exists(os.path.join(outDir,"models")):
    os.makedirs(os.path.join(outDir,"models"))

from pycisTopic.lda_models import *
model=evaluate_models(models,
                     select_model=90,
                     return_model=True,
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + 'models/model_selection.pdf')

# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

# Save
with open(outDir + 'cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)
