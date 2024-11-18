#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=databasebuild.out
#SBATCH --job-name=build1

module load anaconda3
source activate snpeff

cd /storage/home/tlc458/.conda/envs/snpeff/share/snpeff-5.1-2


java -jar  snpEff.jar build -noCheckCds -noCheckProtein -gff3 -c /storage/home/tlc458/.conda/envs/snpeff/share/snpeff-5.1-2/snpEff.config -v Apal 
