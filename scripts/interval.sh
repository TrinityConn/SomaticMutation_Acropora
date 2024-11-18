#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=10gb
#SBATCH --output=interval.out
#SBATCH --job-name=interval

cd ~/scratch

module load picard 

picard ScatterIntervalsByNs REFERENCE= ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta  OUTPUT_TYPE=ACGT OUTPUT=~/scratch/apal.interval_list



