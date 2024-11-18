#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=bwa_5847

##load bwa
module load bwa

module load samtools

module load picard 

##create directory with indexed files

bwa mem -t 10  ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta ~/scratch/5847_1P.fastq.gz ~/scratch/5847_2P.fastq.gz > ~/scratch/5847_pairedends.sam

samtools view -S -b ~/scratch/5847_pairedends.sam > ~/scratch/5847_final.bam

samtools sort  ~/scratch/5847_final.bam -o ~/scratch/5847_sorted_final.bam


picard AddOrReplaceReadGroups I=~/scratch/5847_sorted_final.bam O=~/scratch/5847_vready.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=colony2 RGSM=5847

picard MarkDuplicates I=~/scratch/5847_vready.bam O=~/scratch/5847_dup_final.bam M=~/scratch/marked_dup_5847.txt REMOVE_SEQUENCING_DUPLICATES=true

samtools view -q 20 -f 0x2 -b ~/scratch/5847_dup_final.bam > ~/scratch/5847_SoGv.bam

samtools index  ~/scratch/5847_SoGv.bam

samtools flagstat ~/scratch/5847_SoGv.bam > 5847_SoGv_stat.txt
