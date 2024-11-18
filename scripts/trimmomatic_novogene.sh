#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --job-name=trim_132_2

cd $SLURM_SUBMIT_DIR


#run trimmomatic on paired end reads, using a quality sliding window


cd /storage/group/ibb3/default/tools/trimmomatic-master/bin

java -jar trimmomatic.jar PE -threads 4 -phred33 -basein /storage/group/ibb3/default/SomaticMutationDeep/M132_2/M132_2_CKDN220052804-1A_HF33VDSX5_L1_1.fq.gz -baseout ~/scratch/132_2_trimmed.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36
