#!/bin/bash
#PBS -l nodes=5:ppn=1
#PBS -l walltime=48:00:00
#PBS -A open
#PBS -l mem=4gb
#PBS -j oe 
#PBS -N downsamplesamp


module load gcc/8.3.1 
module load parallel/20200822
module load picard 

cd ~/scratch/filesfordup/callfiles/subset
mkdir downsample

picard DownsampleSam I=16631_dup_final.bam  O=~/scratch/filesfordup/callfiles/subset/downsample/16631_downsample.bam STRATEGY=HighAccuracy P=0.7 
