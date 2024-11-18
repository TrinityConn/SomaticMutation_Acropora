#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=varscan71.out
#SBATCH --job-name=varscan71

module load anaconda3

source activate varscan

cd ~/scratch/SoGV_Reruns/alignments
mkdir var71

varscan somatic 5849_SoGv.bam.pileup 5845_SoGv.bam.pileup var71/var71.basename
