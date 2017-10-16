#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=rondror
#SBATCH --tasks=1 --cpus-per-task=1
##SBATCH --ntasks-per-socket=2 --nodes=1
#SBATCH --mem=16GB

RECEPTOR=$1
STRUCTURE=$2

python score.py $RECEPTOR $STRUCTURE
