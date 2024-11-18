#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -l walltime=05:00:00
#PBS -l pmem=10gb
#PBS -A open
#PBS -j oe
#PBS -N hap_int14scaf

conda activate


module load gcc/8.3.1 
module load parallel/20200822 
cd ~/scratch/filesfordup/callfiles/subset/finalfiles

parallel -j 10 gatk  HaplotypeCaller \
		-R ~/scratch/reference/Apalm_assembly_v3.1_200911.masked.fasta \
		-I  {} \
		-O  ~/scratch/filesfordup/gvcfs/{.}_scaffold14.gvcf \
		-mbq 25 \
		--dont-use-soft-clipped-bases \
		--native-pair-hmm-threads 10 \
		-L  hic_scaffold_14 \
		-ERC GVCF ::: *.bam
