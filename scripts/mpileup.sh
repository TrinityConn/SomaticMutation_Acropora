#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=mpileup4.out
#SBATCH --job-name=mpileu4

module load anaconda3

source activate varscan 

cd ~/scratch/SoGV_Reruns/alignments

parallel -j 20 --tagstring {} "samtools mpileup -f ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta {}>{}.pileup"  ::: *SoGv.bam
