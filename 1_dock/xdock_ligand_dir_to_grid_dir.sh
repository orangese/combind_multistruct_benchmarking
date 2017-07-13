#!/bin/bash
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

GRIDS=$1
LIGANDS=$2
OUT=$3

for grid in $(ls $GRIDS)
do
    sbatch --time=20:00:00 --job-name=$grid -n 1 -p rondror /share/PI/rondror/$USER/combind/1_dock/xdock_ligand_dir_to_single_grid.sh $GRIDS/$grid/$grid.zip $LIGANDS $OUT &
done
