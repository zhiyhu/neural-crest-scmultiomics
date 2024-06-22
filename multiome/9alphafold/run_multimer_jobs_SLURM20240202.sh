#!/bin/bash

#A typical run takes couple of hours but may be much longer
#SBATCH --job-name=array
#SBATCH --time=2-00:00:00

#log files:
#SBATCH -e logs/run_multimer_jobs_%A_%a_err.txt
#SBATCH -o logs/run_multimer_jobs_%A_%a_out.txt

#qos sets priority
#SBATCH --qos=low

#SBATCH -p gpu
#lower end GPUs might be sufficient for pairwise screens: ### SBATCH -C "gpu=2080Ti|gpu=3090"

#Reserve the entire GPU so no-one else slows you down
#SBATCH --gres=gpu:1

#Limit the run to a single node
#SBATCH -N 1

#Adjust this depending on the node
#SBATCH --ntasks=8
#SBATCH --mem=64000

eval "$(conda shell.bash hook)"
conda activate AlphaPulldown

module load cuda/12.2

# https://github.com/tensorflow/tensorflow/issues/61468
CUDNN_PATH=$(dirname $(python -c "import nvidia.cudnn;print(nvidia.cudnn.__file__)"))
export LD_LIBRARY_PATH=$CUDNN_PATH/lib:$CONDA_PREFIX/lib/:$LD_LIBRARY_PATH

nvcc --version

nvidia-smi

MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
export TF_FORCE_UNIFIED_MEMORY='1'
 

run_multimer_jobs.py --mode=pulldown \
    --num_cycle=3 \
    --num_predictions_per_model=1 \
    --output_path=multiome/analysis_newref/alphafold/output/models20240202 \
    --data_dir=/data/share_for_user/alphfold_db \
    --protein_lists=multiome/analysis_newref/alphafold/data/bait20240202_2.txt,multiome/analysis_newref/alphafold/data/targets20240202.txt \
    --monomer_objects_dir=multiome/analysis_newref/alphafold/output/features20240128 \
    --job_index=$SLURM_ARRAY_TASK_ID

