#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=indexvcf
#SBATCH --account=ibb3

source ~/.bashrc
cd ~/scratch/ofav_comparisons/ofav101

conda  activate sambcfenv

for F in *_Mutations.vcf.gz.gz ; do  tabix -f -p vcf ${F} ; done 
