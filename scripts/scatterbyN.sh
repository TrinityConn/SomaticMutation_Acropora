#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=scatterofav.out
#SBATCH --job-name=scatterofav


module load picard
cd ~/scratch

picard ScatterIntervalsByNs \
      REFERENCE=um_ofav_v1_softmask.fa \
      OUTPUT_TYPE=ACGT \
      OUTPUT=output.interval_list
