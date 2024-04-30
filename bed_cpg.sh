#!/bin/bash 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=128gb
#SBATCH --output=bedgff1.out
#SBATCH --job-name=cpg1

module load anaconda3

source activate sambcfenv

cd ~/scratch/SoGV_Reruns/filtered

parallel -j 20 --tagstring {} "bedtools intersect -a {} -b ~/scratch/Apal.gff -header -wa > cpg/{.}_cpg.vcf.gz" ::: *_Mutations.vcf.gz
