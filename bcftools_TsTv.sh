#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=tstv_extract
#SBATCH --account open

module load anaconda3

source activate sambcfenv

cd  ~/scratch/

parallel -j 20 --tagstring {} "bcftools query -H -f '%CHROM\t%POS [\t%AF\t%DP]\n' {} > {.}_3samp.txt" ::: *.vcf


