#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=11
#SBATCH --time=48:00:00
#SBATCH --mem=50gb
#SBATCH --job-name=markdup_20033
#SBATCH --account=open
source ~/.bashrc
module load picard

conda activate sambcfenv

cd ~/scratch

picard AddOrReplaceReadGroups I=/storage/group/dut374/default/trinity/bam/20033.bam O=~/scratch/20033_vready.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=colony3 RGSM=20033

picard MarkDuplicates I=~/scratch/20033_vready.bam O=~/scratch/20033_dup_final.bam M=~/scratch/marked_dup_20033.txt REMOVE_SEQUENCING_DUPLICATES=true

samtools view -q 20 -f 0x2 -b ~/scratch/20033_dup_final.bam > ~/scratch/20033_SoGv.bam

samtools index  ~/scratch/20033_SoGv.bam

samtools flagstat ~/scratch/20033_SoGv.bam > 20033_SoGv_stat.txt
