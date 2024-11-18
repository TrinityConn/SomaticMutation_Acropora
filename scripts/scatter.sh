#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=100gb
#SBATCH --output=scatter.out
#SBATCH --job-name=scatter

cd ~/scratch
source ~/.bashrc
conda activate gatk

gatk SplitIntervals -R ~/scratch/um_ofav_v1_softmask.fa -L output.interval_list --scatter-count 30  -O interval-files 

