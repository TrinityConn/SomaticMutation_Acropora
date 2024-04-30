#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=VAF_Ann_extract
#SBATCH --account open

module load anaconda3

source activate sambcfenv

cd  ~/scratch

parallel -j 20 --tagstring {} "bcftools query -H -f '%CHROM\t%POS [\t%AF\t%DP\t%ANN]\n'  {} > {.}_VAF_ann.txt" ::: *.vcf.gz
