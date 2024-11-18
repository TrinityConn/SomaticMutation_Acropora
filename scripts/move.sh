#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=144:00:00
#SBATCH --mem=100gb
#SBATCH --job-name=move
#SBATCH --account=ibb3
#SBATCH --partition=sla-prio
#SBATCH --output=move.out

source ~/.bashrc

#move to directory 

cd /storage/group/dut374/default/trinity/trimmed

cp -r *  ~/scratch/
