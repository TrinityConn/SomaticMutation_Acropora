#!/bin/bash
#SBATCH --account=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=1gb
#SBATCH --job-name=conc_2
#SBATCH --output=bcf_concatofav2.out

cd ~/scratch/ofav_comparisons/ofav85

source ~/.bashrc

conda activate sambcfenv


bcftools concat --file-list ids.txt  --output ~/scratch/ofav_comparisons/ofav85/ofav85.vcf.gz --output-type z

