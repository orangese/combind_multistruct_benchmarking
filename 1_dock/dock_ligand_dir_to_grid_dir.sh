#!/bin/bash
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

GRIDS=$1
LIGANDS=$2
OUT=$3

for grid in $(ls $GRIDS)
do
    sbatch --time=4:00:00 -n 1 -p rondror /share/PI/rondror/docking_code/1_prep_structures/dock_ligand_dir_to_single_grid.sh $GRIDS/$grid/$grid.zip $LIGANDS $OUT &
done
