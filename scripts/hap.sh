#!/bin/bash
#SBATCH --job-name=5849_hap
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=hp_09_hap.out
#SBATCH --error=hp_09_hap.err

module load anaconda3

source activate gatk

cd ~/scratch/SoGV_Reruns/alignments

gatk HaplotypeCaller -R ~/scratch/Apalm_assembly_v3.1_200911.masked.fasta -I 5849_SoGv.bam -O ~/scratch/SoGV_Reruns/hapcall/5849_BP.g.vcf -ERC BP_RESOLUTION

