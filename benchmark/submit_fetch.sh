#!/bin/bash
#
#SBATCH --job-name=fetch
#SBATCH --output=jobs/%j.out

source ../schrodinger_activate
srun python fetch.py feb2023_benchmark.txt ligands
