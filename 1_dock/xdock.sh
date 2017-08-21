#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --partition=rondror
#SBATCH --tasks=5 --cpus-per-task=1
##SBATCH --ntasks-per-socket=2 --nodes=1
#SBATCH --mem=8GB
#SBATCH --job-name=xdock
#SBATCH --output=out_xdock.out --open-mode=append

ACTION=$1
DATASET=$2

./main.py $ACTION $DATASET
