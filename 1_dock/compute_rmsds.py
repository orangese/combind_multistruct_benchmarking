#!/share/PI/rondror/software/miniconda/bin/python

import sys
import os

sys.path.append(sys.argv[2])
from rmsd import RMSD

DATA = '/scratch/PI/rondror/docking_data/'
dataset = sys.argv[1]

OUTPUT = 'rmsd.csv'

os.chdir(DATA+dataset)
grids_dir = DATA+dataset+'/grids'
structs = [d for d in os.listdir(grids_dir) if os.path.isdir(os.path.join(grids_dir, d))]

def references(ligand):
    return "ligands/{}.mae".format(ligand)

def name(ligand, grid):
    return "{}-to-{}".format(ligand, grid)

def glides(name):
    return "glide/{}/{}_pv.maegz".format(name, name)

rmsds = {}

out = open(OUTPUT, 'w')
for grid in structs:
    grid = grid # look out for capitalization issues
    rmsds[grid] = {}
    for ligand in structs:
        ligand = ligand + '_ligand' # look out for capitalization issues
        print references(ligand), glides(name(ligand, grid))
        try:
            rmsds[grid][ligand] = RMSD.create(references(ligand), glides(name(ligand, grid)))
        except IOError:
            rmsds[grid][ligand] = RMSD([10.0])
        out.write(name(ligand, grid) + ':' + str(rmsds[grid][ligand]) + '\n')

