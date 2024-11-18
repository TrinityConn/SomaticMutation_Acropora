#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem=100gb
#SBATCH --output=gwas.out
#SBATCH --job-name=filter
#SBATCH --account=open

sourcen ~/.bashrc

mamba activate plink2

cd ~/scratch/AGF_GWAS


plink2 --bfile AGF --pheno AGF_Phenotype_Matrix.txt --glmm  --covar plink2.eigenvec  --king-cutoff 0.20 --out agf_gwas  
