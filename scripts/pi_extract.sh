#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=pi_extract
#SBATCH --account open

module load anaconda3

source activate sambcfenv

cd  /storage/group/ibb3/default/HSFP_DeepSeq/final_vcfs

mkdir ~/scratch/pi_window1

parallel -j 20 --tagstring {} "vcftools --gzvcf {} --window-pi 1 --out ~/scratch/pi_window1/{.}" ::: *_Mutations.vcf.gz*
