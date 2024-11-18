#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=48:00
#SBATCH --mem=10gb
#SBATCH --job-name=genome_dict

cd ~/scratch 
module load picard 


picard CreateSequenceDictionary \
	R=um_ofav_v1_softmask.fa
	O=um_ofav_v1_softmask.fa.dict


