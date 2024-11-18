#!/bin/bash
#PBS -l nodes=1:ppn=10:himem
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=100gb
#PBS -j oe
#PBS -N validate_20


module load picard

cd ~/scratch/bamfiles 

picard ValidateSamFile I=16620.bam MODE=SUMMARY
