#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=144:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=ofav1
#SBATCH --account=ibb3
#SBATCH --partition=sla-prio
#SBATCH --output=ofav1.out
#SBATCH --error=ofav1.err

##load bwa
source ~/.bashrc
conda  activate sambamba

##create directory with indexed files

cd ~/scratch

bwa mem -t 20  ~/scratch/combined_ref.fa    ~/scratch/20010_1trimmed.fq.gz ~/scratch/20010_2trimmed.fq.gz |samtools view -S -b  -o - - | samtools sort -o 20010_sym.bam


