#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=10:00:00
#SBATCH --mem=128gb
#SBATCH --output=table.out
#SBATCH --job-name=table
#SBATCH --account open 

source ~/.bashrc

conda activate gatk 

cd /storage/group/ibb3/default/HSFP_DeepSeq/snpeff_results/vcfs


gatk VariantsToTable \
	-V aged92_snpEff.vcf.gz.gz \
	--output aged92.table
