#!/bin/bash

STRUCT=$1
RECEPTOR=$2
SCHRODINGER=/share/PI/rondror/software/schrodinger2017-1

cd $RECEPTOR/processed/

$SCHRODINGER/utilities/prepwizard -WAIT -fillsidechains -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini -j temp-$STRUCT $STRUCT'_before.mae' $STRUCT'_after.mae'
