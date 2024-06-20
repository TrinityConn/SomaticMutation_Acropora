#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=10:00:00
#SBATCH --mem=128gb
#SBATCH --output=C4.out
#SBATCH --job-name=C4
#SBATCH --account ibb3

source ~/.bashrc

conda activate gatk 

mkdir ~/scratch/colony15_4
cd ~/scratch

parallel -j 20 gatk Mutect2 -R Apalm_assembly_v3.1_200911.masked.fasta \
 -I ~/scratch/S15_1_SoGv.bam \
 -I ~/scratch/S15_5_SoGv.bam \
 -tumor S15_5 \
 -normal S15_1 \
 -O ~/scratch/colony15_4/c15_4_{.}.vcf.gz \
 -L {} ::: *.intervals
