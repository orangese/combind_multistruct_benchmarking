#!/share/PI/rondror/software/miniconda/bin/python

import sys
import os

sys.path.append(sys.argv[2])
from rmsd import RMSD

DATA = '/scratch/PI/rondror/docking_data/'
dataset = sys.argv[1]

os.chdir(DATA+dataset)

def references(ligand):
    return "ligands/{}_ligand.mae".format(ligand)

def name(ligand, grid):
    return "{}_ligand-to-{}".format(ligand, grid)

def xglides(name):
    return "xglide/{}/{}_pv.maegz".format(name, name)

def glides(name):
    return "glide/{}/{}_pv.maegz".format(name, name)

for i, dock in enumerate(['/glide/']):#, '/xglide/']):
    rmsds = {}

    gdir = DATA+dataset+dock
    if not os.path.isdir(gdir): continue

    if i == 0: out = open('rmsd.csv', 'w')
    if i == 1: out = open('xrmsd.csv', 'w')

    for fnm in os.listdir(gdir):
        lig, struct = fnm.split('_ligand-to-')
        
        if not (os.path.exists(gdir + fnm + '/' + fnm + '_pv.maegz') and os.path.exists(gdir + fnm + '/' + fnm + '.rept')):
            continue
        
        if struct not in rmsds: rmsds[struct] = {}
        
        if i == 0: loc = glides(name(lig, struct))
        if i == 1: loc = xglides(name(lig, struct))

        try:
            rmsds[struct][lig] = RMSD.create(references(lig), loc)
        except IOError:
            rmsds[struct][lig] = RMSD([10.0])
        out.write(name(lig, struct) + ':' + str(rmsds[struct][lig]) + '\n')

