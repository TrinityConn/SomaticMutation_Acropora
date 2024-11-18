
#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=10gb
#PBS -j oe
#PBS -N genotypevcf


cd ~/scratch/filesfordup

conda activate

gatk GenotypeGVCFs -R /gpfs/group/ibb3/default/HSFP_DeepSeq/alignments/apal.fasta -V 16629.g.vcf.gz -O final_youngtest.vcf.gz

