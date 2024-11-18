#!/bin/bash
#SBATCH --job-name=emboss_cpg
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBAATCH --output=cpg_apal.out
#SBATCH --error=cpg_apal.err

cd ~/scratch

module load anaconda3

source activate gatk 

cpgreport -sequence Apalm_assembly_v3.1_200911.masked.fasta -score 28 -outfile Apal_cpg -outfeat Apal.gff
