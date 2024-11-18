#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:0:00
#SBATCH --mem=100gb
#SBATCH --output=gffread.out
#SBATCH --job-name=gffread


module load anaconda3

source activate gffread

cd ~/scratch

gffread -g Apalm_assembly_v3.1_200911.masked.fasta genes.gff -x Apal_CDS.fasta
