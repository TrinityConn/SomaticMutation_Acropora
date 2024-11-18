#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=20gb
#PBS -j oe 
#PBS -N mutect21opp


cd ~/scratch

gatk Mutect2 -R ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta \
 -I dup_16628.bam \
 -I dup_16627.bam \
 -normal 16627 \
 -O somatic_2.vcf.gz \
 -L ~/work/interval-files-folder/0000-scattered.interval_list
