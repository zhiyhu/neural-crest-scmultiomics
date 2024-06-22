#!/bin/sh
#SBATCH --job-name=ctx_bd
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log 
## Job needs 1 nodes and 4 cores per node.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # specify number of CPUs to use here
## Job request memory
#SBATCH --mem-per-cpu=15gb

## Zhiyuan Hu
## 19 Nov 2022
## last modified 14 Dec 2023

## multiome/analysis_newref/GRN_scenicplus/ncall_2023oct_ccb/code/3create_cisTarget_databases.sh

## Prior condition:
### 1.Run Project_analysis/multiome/analysis_newref/GRN_scenicplus/ncall/code/prepare_motif2tf.R
### 2.Run /Filers/home/z/zhu/t1data/multiome/analysis_newref/GRN_scenicplus/motif_process/code/prepare_cbuster_motifs_newm.py

proj=ncall_2023oct_ccb
ncpu=4

wkdir=multiome/analysis_newref/GRN_scenicplus
# Paths and parameters
consensdir=${wkdir}/${proj}/output/consensus_peak_calling
outdir=${wkdir}/${proj}/output/ctx_db # cisTopic_database
tag='ncall_newmotif' ## tag
genomefa='ref/ensembl105/GRCz11.105_foxd3_mcherry_citrine/Danio_rerio.GRCz11_foxd3_mcherry_citrine.dna.primary_assembly.fa'

cbdir=${wkdir}/data/motif_with_newm/motifs_cb_format #Path to directory with Cluster-Buster motifs
# Create outdir
if [ ! -d $outdir ]
then
    mkdir -p $outdir
fi

#### List of motifs
motif_list=${wkdir}/data/motif_with_newm/motifs.lst

#### Get fasta sequences
echo "Extracting FASTA ..."
module load bedtools/2.29.2
bedtools getfasta -fi $genomefa \
                  -bed $consensdir/consensus_regions.bed > \
                  $consensdir/consensus_regions.fa
echo "Done."

#### Create scores DB
echo "Creating scores DB files ..."
#### Activate environment
module purge
source /home/z/zhu/miniconda3/etc/profile.d/conda.sh
source /home/z/zhu/miniconda3/bin/activate
conda activate .conda/envs/create_cistarget_databases

conda info
conda list flatbuffers
#### Set ${create_cistarget_databases_dir} to https://github.com/aertslab/create_cisTarget_databases
create_cistarget_databases_dir='multiome/R/GRN_SCENICplus/create_cisTarget_databases'

#### Score the motifs in 1 chunks; we will use the non-redundant db here
# for current_part in {1..10} ; do
.conda/envs/create_cistarget_databases/bin/python  ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
         -f $consensdir/consensus_regions.fa \
         -M $cbdir \
         -m $motif_list \
         -o $outdir/$tag \
         -t 35 \
         -l
         # -p ${current_part} 10 \
# done
echo "Done."
#### Create rankings
echo "Creating rankings DB files ..."
.conda/envs/create_cistarget_databases/bin/python  ${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
         -i $outdir/$tag.motifs_vs_regions.scores.feather -s 555
echo "Done."
echo "ALL DONE."

