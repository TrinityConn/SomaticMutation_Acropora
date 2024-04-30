#!/bin/bash
#SBATCH --account=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --out=filterremove.out
#SBATCH --job-name=filter_mutect_parallel

source ~/.bashrc

conda  activate sambcfenv

cd ~/scratch/ofav_comparisons/ofav85


for filename in *_Filtered.vcf.gz; do 
bcftools view -i 'FILTER=="PASS"'  $filename  >   ${filename}_filtered.vcf.gz

done
