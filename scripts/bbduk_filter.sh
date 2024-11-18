#!/bin/bash
#SBATCH --job-name=bbduk20011
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=trim11.out
#SBATCH --account=ibb3

source ~/.bashrc
conda  activate bbmap_trim

cd  /storage/group/dut374/default/trinity/rawdata



bbduk.sh \
in1=/storage/group/dut374/default/trinity/rawdata/Unknown_BK461-020033_1.fq.gz \
in2=/storage/group/dut374/default/trinity/rawdata/Unknown_BK461-020033_2.fq.gz \
out1=~/scratch/20033_1trimmed.fq.gz \
out2=~/scratch/20033_2trimmed.fq.gz \
ref=/storage/group/ibb3/default/tools/bbmap/resources/adapters.fa \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
tpe \
tbo \
qtrim=r \
tossbrokenreads=t \
trimq=20 \
maq=20 \
minlen=50 
