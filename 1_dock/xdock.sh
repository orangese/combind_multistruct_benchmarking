#!/bin/sh

#SBATCH --time=10:00:00
#SBATCH --partition=rondror

ACTION=$1
DATASET=$2

./main.py $ACTION $DATASET
