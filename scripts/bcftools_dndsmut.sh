#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=vaf_extract
#SBATCH --account open

source ~/.bashrc

conda activate sambcfenv

cd  /storage/group/ibb3/default/HSFP_DeepSeq/final_vcfs

parallel -j 20 --tagstring {} "bcftools query -H -f '%CHROM\t%POS %REF %ALT\n' {} > {.}_dnds.txt" ::: *_Mutations.vcf.gz.gz


