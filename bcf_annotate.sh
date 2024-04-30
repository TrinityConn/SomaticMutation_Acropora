#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=128gb
#SBATCH --output=bcfannotate.out
#SBATCH --job-name=annotate1

module load anaconda3
source activate sambcfenv
mkdir ~/scratch/aged/final_vcfs/filtered_prelim/final/bed_out

cd ~/scratch/aged/final_vcfs/filtered_prelim/final

for filename in *_Mutations.vcf.gz; do 
 vcftools --gzvcf ${filename} --bed ~/scratch/Apalm.bed --out bed_out/${filename}_bed --recode 
done

