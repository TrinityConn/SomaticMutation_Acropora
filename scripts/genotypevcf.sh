#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=10gb
#PBS -j oe
#PBS -N genotypevcf


cd ~/scratch/filesfordup

conda activate 


gatk CombineGVCFs -R /gpfs/group/ibb3/default/HSFP_DeepSeq/alignments/apal.fasta --variant 16620.g.vcf.gz --variant 16622.g.vcf.gz  --variant 16627.g.vcf.gz --variant 16625.g.vcf.gz --variant 16628.g.vcf.gz --variant 16629.g.vcf.gz --o young colony.g.vcf.gz


gatk GenotypeGVCFs -R /gpfs/group/ibb3/default/HSFP_DeepSeq/alignments/apal.fasta -V youngcolony.g.vcf.gz -O final_young.vcf.gz 




