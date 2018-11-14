#!/bin/bash
#SBATCH --partition=rondror
#SBATCH --time=02:00:00
#SBATCH --open-mode=truncate
#SBATCH --mem=4GB
#SBATCH --out=/home/users/jpaggi/combind/tests/outputs/glide_ligands.out
#SBATCH --err=/home/users/jpaggi/combind/tests/outputs/glide_ligands.out

. setup_combind

$SCHRODINGER/run $COMBINDHOME/dock/prepare_all.py 'g'
