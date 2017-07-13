#!/bin/bash
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GSCRIPT=/share/PI/rondror/$USER/combind/1_dock/xglide_SP.sh
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
