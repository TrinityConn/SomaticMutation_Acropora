#!/bin/bash
#SBATCH --job-name=bbduk16628
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=16629_repair.out
#SBATCH --account=ibb3_d

module load anaconda3

source activate bbmap_trim

cd /storage/group/ibb3/default/SomaticMutationDeep/young_colony/

repair.sh \
in1=16629_S3_R1_001.fastq.gz \
in2=16629_S3_R2_001.fastq.gz \
out1=16629_R1_fixed.fastq.gz \
out2=16629_R2_fixed.fastq.gz \
outs=singletons.fq repair

