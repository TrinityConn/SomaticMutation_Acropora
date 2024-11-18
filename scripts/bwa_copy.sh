#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=144:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=20033
#SBATCH --account=ibb3
#SBATCH --partition=sla-prio
#SBATCH --output=20033.out
#SBATCH --error=20033.err

##load bwa
source ~/.bashrc
conda  activate sambamba

##create directory with indexed files

cd ~/scratch

bwa mem -t 20  ~/scratch/ofav_comparisons/um_ofav_v1_softmask.fa   ~/scratch/20033_1trimmed.fq.gz ~/scratch/20033_2trimmed.fq.gz |samtools view -S -b  -o - - | samtools sort -o 20033.bam



