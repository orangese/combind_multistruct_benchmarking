#!/bin/bash
##SBATCH --time=4:00:00
##SBATCH --partition=rondror
#SBATCH --cpus-per-task=1
##SBATCH --tasks=2 --cpus-per-task=1
#SBATCH --ntasks-per-socket=2
#SBATCH --mem=16GB

PROTEIN=$1
QUERY=$2
python score_script.py $PROTEIN $QUERY $3 $4 $5
