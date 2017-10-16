#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --partition=rondror
#SBATCH --tasks=2 --cpus-per-task=1
##SBATCH --ntasks-per-socket=2 --nodes=1
#SBATCH --mem=16GB

#module load chemistry schrodinger/2017-2

ACTION=$1
DATASET=$2

$SCHRODINGER/run main.py $ACTION $DATASET
