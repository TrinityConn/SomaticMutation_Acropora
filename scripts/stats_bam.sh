#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=10gb
#SBATCH --output=pileup_20033.out
#SBATCH --job-name=pileup20033
#SBATCH --account ibb3

source ~/.bashrc
conda activate bbmap_trim
pileup.sh in=/storage/group/dut374/default/trinity/bam/20033_SoGv.bam  ref=~/scratch/ofav_comparisons/um_ofav_v1_softmask.fa  out=~/scratch/20033
