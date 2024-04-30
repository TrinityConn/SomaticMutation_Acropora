#!/bin/bash
#SBATCH --account=open
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=filterofav2

source ~/.bashrc

module load samtools

conda activate gatk 

cd ~/scratch/ofav_comparisons/ofav85

parallel -j 20 gatk FilterMutectCalls \
	-V  {} \
	 -O {.}_Filtered.vcf.gz \
	 -R ~/scratch/ofav_comparisons/um_ofav_v1_softmask.fa \
	--min-median-base-quality 25 \
	--min-median-mapping-quality 30 ::: ofav85.vcf.gz 



