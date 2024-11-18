#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=runmanta1.out
#SBATCH --job-name=manta_run


module load anaconda3

source activate sambcfenv 

cd ~/scratch/SoGV_Reruns/alignments/manta2

./runWorkflow.py -j 8 
