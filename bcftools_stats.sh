#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=1:00:00
#SBATCH --mem=10gb
#SBATCH --out=stats.out
#SBATCH --job-name=extract_stats
#SBATCH --account=ibb3

source ~/.bashrc
conda activate sambcfenv



cd  ~/scratch

parallel -j 20 --tagstring {} "bcftools stats {} > {.}_stats.txt" ::: ofav*_Mutations.vcf.gz.gz






