#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --job-name=132_2
#SBATCH --account=open
#SBATCH --output=bwa_132_2.out

##load bwa
module load bwa

module load samtools

module load picard 

##create directory with indexed files
cd ~/scratch
bwa index um_ofav_v1_hardmask.fa

bwa mem -t 10  ~/scratch/um_ofav_v1_hardmask.fa  ~/scratch/10010_1trimmed.fq.gz ~/scratch/10010_2trimmed.fq.gz > ~/scratch/10010_pairedends.sam
