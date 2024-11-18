#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=trim_5849

cd $SLURM_SUBMIT_DIR


#run trimmomatic on paired end reads, using a quality sliding window


cd /storage/group/ibb3/default/tools/trimmomatic-master/bin
java -jar trimmomatic.jar  PE -threads 4  -phred33  -basein /storage/group/ibb3/default/HSFP_DeepSeq/rawdata/old/5849_S20_R1_001.fastq.gz   -baseout ~/scratch/5849_trimmed.fastq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
