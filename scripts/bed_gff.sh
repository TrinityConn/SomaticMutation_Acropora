#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=bedgenes.out
#SBATCH --job-name=gene1

module load anaconda3

source activate sambcfenv

cd ~/scratch/SoGV_Reruns/filtered

parallel -j 20 --tagstring {} "bedtools intersect -a {} -b /storage/group/ibb3/default/AP_AC_genome_seqs/dovetail_Apalm/HiC_improvement/Apalm_assembly_v3.1_200911.gff3  -header -wa -wb  > {.}_genes.vcf.gz" ::: *_Mutations.vcf.gz
