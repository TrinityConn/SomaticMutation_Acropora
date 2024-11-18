#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --job-name=extract_st
#SBATCH --out=extract_st.out
#SBATCH --time=00:30:00
#SBATCH --account=ibb3



source ~/.bashrc
conda activate sambcfenv

cd ~/scratch/

parallel -j 20 --tagstring {} "awk '$1 /^ST/' {} > {.}_ST.txt" ::: ofav*_stats*


