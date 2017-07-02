#!/bin/bash

RECEPTOR=$1
SCHRODINGER=/share/PI/rondror/software/schrodinger2017-1
DATA=/scratch/PI/rondror/docking_data/

SCRIPT=/share/PI/rondror/$USER/combind/1_dock/schro_process.sh

if [ ! -d $DATA$RECEPTOR/processed ]
then
    mkdir $DATA$RECEPTOR/processed
fi

for struct in $(ls $DATA$RECEPTOR/stripped)
do
    if [ ! -e $DATA$RECEPTOR/processed/$struct ]
    then
        struct=${struct:0:4}
        cp $DATA$RECEPTOR/stripped/$struct'.mae' $DATA$RECEPTOR/processed/$struct'_before.mae'
        echo Processing $struct...
        sbatch --time=2:00:00 --output=$DATA$RECEPTOR/processed/slurm-%A_%a.out --job-name=p-$RECEPTOR -n 1 -p rondror $SCRIPT $struct $DATA$RECEPTOR
    fi
done
