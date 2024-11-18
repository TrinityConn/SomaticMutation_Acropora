#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --job-name=filter

module load anaconda3

source activate gatk 
cd ~/scratch

gatk VariantFiltration \
-R Apalm_assembly_v3.1_200911.masked.fasta \
-V Old1_Final.vcf.gz \
-O Old1_filter1.vcf.gz \ 
-select 'vc.getGenotype("5838").getDP() > 12'
