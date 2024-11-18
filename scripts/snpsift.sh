#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=10:00:00
#SBATCH --mem=10gb
#SBATCH --out=snpsift.out

source ~/.bashrc



conda activate snpeff


cd /storage/group/ibb3/default/HSFP_DeepSeq/snpeff_results/vcfs

for filename in *.gz.gz; do

SnpSift filter "EFF[*].EFFECT='SYNONYMOUS' ${filename}"  > ${filename}_syn.vcf
done

