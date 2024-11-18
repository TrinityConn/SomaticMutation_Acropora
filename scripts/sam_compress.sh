#!/bin/bash
#PBS -l nodes=1:ppn=10:himem
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N comp_16629


cd ~/scratch

samtools view -S -b 16629_pairedends_final.sam > 16629_final.bam
