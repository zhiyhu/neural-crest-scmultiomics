#!/bin/bash
#SBATCH --job-name=prepbam
#SBATCH -o logs/%j_%x.log 
#SBATCH -e logs/%j_%x.log
#SBATCH --nodes=1
#SBATCH --mem=30gb

## Zhiyuan
## 18 OCT 2023
## last modified 18 OCT 2023

eval "$(conda shell.bash hook)"

conda activate chrombpnet

cluster=$1

echo $cluster

bam=chrombpnet/data/03merged_bam/${cluster}.bam
new_bam=chrombpnet/data/03merged_bam/${cluster}_new.bam

module load samtools

mkdir chrombpnet/data/tmp/${cluster}
cd chrombpnet/data/tmp/${cluster}

# Convert BAM to SAM
samtools view -H $bam > header.sam
samtools view $bam > body.sam

# Edit the header section
awk -F '\t' -v OFS='\t' '/^@SQ/ {sub(/^SN:/, "SN:chr", $2)} 1' header.sam > edited_header.sam

# Edit the alignment section
awk -F '\t' -v OFS='\t' '{$3="chr"$3; if ($7 ~ /^[0-9XYM]+$/) $7="chr"$7} 1' body.sam > edited_body.sam

# Combine the edited header and body sections
cat edited_header.sam edited_body.sam > edited.sam

# Convert edited SAM back to BAM
# samtools view -b -o $newbam edited.sam
samtools view -h -b -S edited.sam > $new_bam

# Clean up intermediate files (optional)
rm header.sam body.sam edited_header.sam edited_body.sam edited.sam
