#!/bin/bash
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GSCRIPT=/scratch/PI/rondror/docking-julia/1_prep_structures/glide_SP.sh
grid=$1
LIGANDS=$2
OUT=$3

for ligand in $(ls $LIGANDS)
do
    l=${ligand%.*}
    g=${grid%.*}
    l=${l##*/}
    g=${g##*/}
    sh $GSCRIPT $grid $LIGANDS/$ligand $l-to-$g $OUT
done
