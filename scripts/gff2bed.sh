#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --time=2:00:00
#SBATCH --mem=100gb
#SBATCH --output=gff2bed.out
#SBATCH --job-name=gff2bed


module load anaconda3
source activate sambcfenv


cd ~/scratch/SoGV_Reruns/filtered


bedmap --echo --echo-map-id-uniq <(vcf2bed < Old1_Final.vcf_Mutations.vcf.gz) <(gff2bed < ~/scratch/Apalm_assembly_v3.1_200911.gff3) > old1_annotate.bed

