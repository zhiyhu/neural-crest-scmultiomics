#!/bin/bash
#SBATCH --job-name=tfmodisco
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80gb
#SBATCH --time=14-00:00:00

## Zhiyuan
## 18 OCT 2023
## last modified 12 Dec 2023

module load cuda/11.2
eval "$(conda shell.bash hook)"
conda activate /home/huzhiy/miniforge3/envs/tfmodiscolit

cluster=$1  #  mNC_nohox #  # mNC_nohox # Mutant_nohox_12_22ss
echo $cluster
wkdir=/home/huzhiy/projects_ox/multiome/analysis_newref/chrombpnet
outdir=${wkdir}/data/10tfmodisco/${cluster}
h5py=${wkdir}/data/09contribs_bw/${cluster}.profile_scores.h5
mkdir -p $outdir
modisco motifs -i ${h5py} -w 500 -n 1000000 -v -o ${outdir}/modisco_results_profile.h5

# contribution score profiles derived from a ChromBPNet model
# cropped to 500 bp around the peak summit
# examples: https://github.com/kundajelab/chrombpnet/blob/8e7fd8ae71eaadae6594a98c18ee7c93ebd6eb2e/workflows/train_chrombpnet_model.sh#L428

# options:
#   -h, --help            show this help message and exit
#   -s SEQUENCES, --sequences SEQUENCES
#                         A .npy or .npz file containing the one-hot encoded sequences.
#   -a ATTRIBUTIONS, --attributions ATTRIBUTIONS
#                         A .npy or .npz file containing the hypothetical attributions, i.e., the attributions for all nucleotides at all positions.
#   -i H5PY, --h5py H5PY  A legacy h5py file containing the one-hot encoded sequences and shap scores.
#   -n MAX_SEQLETS, --max_seqlets MAX_SEQLETS
#                         The maximum number of seqlets per metacluster.
#   -l N_LEIDEN, --n_leiden N_LEIDEN
#                         The number of Leiden clusterings to perform with different random seeds.
#   -w WINDOW, --window WINDOW
#                         The window surrounding the peak center that will be considered for motif discovery.
#   -o OUTPUT, --output OUTPUT
#                         The path to the output file.
#   -v, --verbose         Controls the amount of output from the code.
