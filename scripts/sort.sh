#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=sort_38

module load samtools

samtools sort  ~/scratch/5838_final.bam -o ~/scratch/5838_sorted_final.bam


module load picard

picard AddOrReplaceReadGroups I=~/scratch/5838_sorted_final.bam O=~/scratch/5838_vready.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=colony1 RGSM=5838

picard MarkDuplicates I=~/scratch/5838_vready.bam O=~/scratch/5838_dup_final.bam M=~/scratch/marked_dup_5838.txt REMOVE_SEQUENCING_DUPLICATES=true

