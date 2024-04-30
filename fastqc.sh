#!/bin/bash 
#PBS -l nodes=1:ppn=10
#PBS -l walltime=48:00:00
#PBS -A open 
#PBS -l mem=10gb
#PBS -j oe
#PBS -N fastqc_16625t

cd $PBS_O_WORKDIR
module load fastqc
module load gcc/8.3.1 
module load parallel/20200822 


cd ~/scratch/trimmed

fastqc 16625_trimmed_1P.fastq.gz   


echo 'done'
