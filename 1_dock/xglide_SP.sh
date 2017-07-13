#!/bin/bash

# Example usage:
# sh glide_SP.sh 3eml-basic-grid/3eml-basic-grid.zip ZMA.mae ZMA_to_3eml

GRID=$1
LIGAND=$2
OUT=$3
GDIR=$4
SCHRODINGER=/share/PI/rondror/software/schrodinger2017-1
cd /share/PI/rondror/$USER/combind/1_dock/
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo $OUT

if [ -e $GDIR$OUT/$OUT'_pv.maegz' ]
then
    echo $OUT already successfully docked
else
    if [ -d $GDIR$OUT ]
    then
        echo deleting failed docking attempt for $OUT
        rm -r $GDIR$OUT
    fi

    echo attempting to dock $OUT

    mkdir $GDIR$OUT
    cp $GRID $GDIR$OUT/$OUT'_grid.zip'
    cp $LIGAND $GDIR$OUT/$OUT'_ligand.mae'

    ################################################################################
    cd $GDIR$OUT

    # Create Glide input file
    sh /share/PI/rondror/$USER/combind/1_dock/xSP-in-file.sh $OUT > $OUT.in
    # Run Glide
    $SCHRODINGER/glide $OUT.in -WAIT

    cd -
fi
