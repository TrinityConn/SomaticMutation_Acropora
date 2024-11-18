#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=index1.out
#SBATCH --job-name=indexout

module load anaconda3

source activate sambcfenv 

cd /storage/group/dut374/default/trinity/bam/20026

for filename in *.bam; do 

samtools index ${filename}


done
