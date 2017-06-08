#!/bin/bash

# Example usage:
# sh glide_SP.sh 3eml-basic-grid/3eml-basic-grid.zip ZMA.mae ZMA_to_3eml

GRID=$1
LIGAND=$2
OUT=$3
GDIR=$4
SCHRODINGER=/share/PI/rondror/software/schrodinger2016-1
cd /share/PI/rondror/docking_code/1_dock/
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo Writing docking results for ligand $LIGAND to grid $GRID to directory $GDIR$OUT

mkdir $GDIR$OUT
cp $GRID $GDIR$OUT/$OUT'_grid.zip'
cp $LIGAND $GDIR$OUT/$OUT'_ligand.mae'

################################################################################
cd $GDIR$OUT

# Create Glide input file
sh /share/PI/rondror/docking_code/1_dock/SP-in-file.sh $OUT > $OUT.in
# Run Glide
$SCHRODINGER/glide $OUT.in -WAIT

cd -
