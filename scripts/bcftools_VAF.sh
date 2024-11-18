#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=vaf_extract
#SBATCH --account open

module load anaconda3

source activate sambcfenv

cd  ~/scratch/ofav_comparisons/ofav86

parallel -j 20 --tagstring {} "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT [\t%AF\t%DP]\n' {} > {.}_refaltsamp.txt" ::: *_Mutations.vcf.gz.gz
