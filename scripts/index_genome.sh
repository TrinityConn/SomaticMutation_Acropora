#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=index_genome
#SBATCH --out=index.out

source ~/.bashrc

conda activate sambcfenv 


cd ~/scratch


samtools faidx clade_ABCDFG_concat.fasta.gz.gz
samtools dict clade_ABCDFG_concat.fasta.gz.gz
