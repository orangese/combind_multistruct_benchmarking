#!/bin/bash
#
#SBATCH --job-name=filter
#SBATCH --output=jobs/%j.out

source ../schrodinger_activate
srun python filter.py ligands 20 additional_smiles
