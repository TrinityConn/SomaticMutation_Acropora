#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=manta276.out
#SBATCH --job-name=manta276

module load anaconda3
source activate sambcfenv

cd ~/scratch/aged/alignments

configManta.py  --normalBam /storage/group/ibb3/default/HSFP_DeepSeq/M132_5_SoGv.bam \
 --tumorBam /storage/group/ibb3/default/HSFP_DeepSeq/M132_4_SoGv.bam \
 --referenceFasta ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta \
 --runDir ./manta276

cd ~/scratch/aged/alignments/manta276

./runWorkflow.py 
