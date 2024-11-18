#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=10:00:00
#SBATCH --mem=128gb
#SBATCH --output=snpeff_young56.out
#SBATCH --job-name=annotagedyoung56
module load anaconda3
source activate snpeff


cd ~/scratch/aged/final_vcfs/filtered_prelim/final


snpEff -v  Apal \
young56_Final.vcf_Mutations.vcf.gz.gz \
-stats young56 \
> young56_snpEff.vcf.gz

