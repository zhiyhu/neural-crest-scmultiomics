#!/bin/sh
#SBATCH --job-name=topic_mod

## Zhiyuan Hu
## 17 Nov 2022
## last modified 31 Dec 2022

## refer to: https://pycistopic.readthedocs.io/en/latest/Cortex_pycisTopic.html#5.-Run-models

i=$1
proj=$2

echo "Run Topic Modelling with nt:"
echo $i


# module load java/17.0.1 
module load mallet

tmp=/t1-data/project/tsslab/zhu/tmp/runModels/nt${i}
wkdir=/t1-data/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/${proj}/output
outdir=${wkdir}/models
out=${outdir}/mallet_nt${i}.pkl
in=${wkdir}/cisTopicObject.pkl

if [ ! -d "$wkdir" ]; then
  mkdir $wkdir
fi

if [ ! -d "$outdir" ]; then
  mkdir $outdir
fi

if [ ! -d "$tmp" ]; then
  mkdir $tmp
fi

cd /t1-data/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/ncall/
python /t1-data/project/tsslab/zhu/multiome/analysis_newref/GRN_scenicplus/nc5k/code/runModels_lda_mallet.py \
        -i $in \
        -o $out \
        -nt $1 \
        -c 4 \
        -it 500 \
        -a 50 \
        -abt True \
        -e 0.1 \
        -ebt False \
        -s 555 \
        -sp $tmp \
        -td $tmp
        
