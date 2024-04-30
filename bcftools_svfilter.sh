#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=vaf_extract
#SBATCH --account open

module load anaconda3

source activate sambcfenv

cd  ~/scratch/aged/alignments/mantaresults

parallel -j 20 --tagstring {} "bcftools query -H -f '%CHROM\t%POS [\t%SVTYPE\t%SVLEN\t%CIPOS\t%CIEND\t%MATEID\t%EVENT\t%SVINSLEN\t%SOMATICSCORE]\n' {} > {.}_SVinfo.txt" ::: *_filtered.vcf.gz


