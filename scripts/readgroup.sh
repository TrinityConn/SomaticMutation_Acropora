#!/bin/bash
#PBS -l nodes=1:ppn=10:himem
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N rg_16625


cd ~/scratch/filesfordup

module load picard 

picard AddOrReplaceReadGroups I=16625_sorted.bam O=16625_vready.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=colony1 RGSM=16625
