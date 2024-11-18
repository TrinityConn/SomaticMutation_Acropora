#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=10010
#SBATCH --account=open
#SBATCH --output=10010.out

##load bwa
module load bwa

module load samtools

module load picard 

##create directory with indexed files

cd ~/scratch 

bwa mem -t 16  ~/scratch/um_ofav_v1_softmask.fa  ~/scratch/10010_1trimmed.fq.gz ~/scratch/10010_2trimmed.fq.gz >  ~/scratch/10010_pairedends.sam




