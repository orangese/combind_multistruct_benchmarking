#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --partition=rondror
#SBATCH --tasks=4 --cpus-per-task=1
##SBATCH --ntasks-per-socket=2 --nodes=1
#SBATCH --mem=16GB
#SBATCH --job-name=xdock
##SBATCH --output=$2.out --open-mode=append

module load chemistry schrodinger/2017-2

ACTION=$1
DATASET=$2

./main.py $ACTION $DATASET
