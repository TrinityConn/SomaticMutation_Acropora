#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l pmem=10gb
#PBS -A open
#PBS -j oe
#PBS -N bcf_merge

conda activate

cd /gpfs/group/ibb3/default/HSFP_DeepSeq/

bcftools merge Somatic_18_filteredfinal.vcf.gz Somatic_19_filteredfinal.vcf.gz Somatic_31_filteredfinal.vcf.gz Somatic_39_filteredfinal.vcf.gz Somatic_43_filteredfinal.vcf.gz Somatic_45_filteredfinal.vcf.gz Somatic_46_filteredfinal.vcf.gz --force-samples -m none > All.vcf

