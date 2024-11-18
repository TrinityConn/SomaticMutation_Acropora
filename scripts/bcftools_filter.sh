#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=20:00:00
#SBATCH --mem=10gb
#SBATCH --output=filter.out
#SBATCH --job-name=filter
#SBATCH --account=open

source ~/.bashrc
conda activate sambcfenv

cd ~/scratch/ofav_comparisons/ofav85


parallel -j 20 --tagstring {} "bcftools view --max-alleles 2 --exclude-types indels,mnps  {} > {.}_Mutations.vcf.gz" ::: ofav85.vcf_Filtered.vcf.gz_filtered.vcf.gz
