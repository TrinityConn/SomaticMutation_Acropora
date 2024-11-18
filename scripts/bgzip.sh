#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=1:00:00
#SBATCH --mem=12gb
#SBATCH --output=bgzip.out
#SBATCH --job-name=bgzip
#SBATCH --account=ibb3

source ~/.bashrc

conda activate sambcfenv

cd  ~/scratch/ofav_comparisons/ofav101

for filename in *_Mutations.vcf.gz; do 

bgzip ${filename} 

done
