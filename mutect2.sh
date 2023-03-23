#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=Lold8.out
#SBATCH --job-name=mutectL8

module load anaconda3 

source activate gatk 

cd ~/scratch

parallel -j 20 gatk Mutect2 -R Apalm_assembly_v3.1_200911.masked.fasta \
 -I ~/scratch/SoGV_Reruns/alignments/5838_SoGv.bam \
 -I ~/scratch/SoGV_Reruns/alignments/5849_SoGv.bam \
 -tumor 5849 \
 -normal 5838 \
 -O ~/scratch/SoGV_Reruns/output/oldL8/oldL8_{.}.vcf.gz \
 -L {} ::: *.interval_list
