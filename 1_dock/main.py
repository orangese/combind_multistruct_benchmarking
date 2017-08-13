#!/share/PI/rondror/software/schrodinger2017-1/run

import sys
import os

import get_pdbs
import stripstructures
import ligprep
import gridgen
import processing
import dock

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
DATA = "/scratch/PI/rondror/docking_data"

commands = sys.argv[1]
datasets = sys.argv[2:]

toRun = []
if 'a' in commands:
    toRun = ['r', 's', 'p', 'l', 'g', 'x']
else:
    toRun = list(commands[1:])

for dataset in datasets:#Go through the given structures, performing the commands specified

    os.chdir('{}/{}'.format(DATA, dataset))

    if 'r' in toRun:
        get_pdbs.get()
    if 's' in toRun:
	stripstructures.strip()
    if 'p' in toRun:
        processing.process()
    if 'l' in toRun:
        ligprep.extractLigands()
    if 'g' in toRun:
        gridgen.getGrids(dataset)
    if 'x' in toRun:
        dock.dockDataset()
