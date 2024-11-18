#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=100gb
#SBATCH --output=pyclone.out
#SBATCH --job-name=pyclone
#SBATCH --account=ibb3


source ~/.bashrc

conda activate pyclone

cd ~/scratch

PyClone run_analysis_pipeline --in_files M1455.tsv M1454.tsv M1453.tsv M1452.tsv M1451.tsv --working_dir ~/scratch/PyClone20


