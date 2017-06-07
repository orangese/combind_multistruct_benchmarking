from rmsd import RMSD
import sys
import os

#structs = ['4ami', '4amj', '3zpr', '3zpq', '2vt4', '2y00', '2y01', '2y02', '2y03', '2y04', '2ycw', '5a8e']
structs = ['4n6h', '4rwa']
#structs = ['4djh', '4dkl', '4ea3', '4n6h', '4rwa', '5c1m', '5dhg', '5dhh']
#structs = ['3rey', '3rfm', '3pwh', '3uza', '3uzc', '2ydv', '2ydo', '4ug2', '3qak']
#structs = ['4dkl', '5c1m']

def references(ligand):
    return "ligands/{}.mae".format(ligand)

def name(ligand, grid):
    return "{}-to-{}".format(ligand, grid)

def glides(name):
    # glide/2VT4_lig-to-2vt4_grid/2VT4_lig-to-2vt4_grid_pv.maegz
    return "glide/{}/{}_pv.maegz".format(name, name)
    
rmsds = {}
if sys.argv[1] == '-o':
    out = open(sys.argv[2], 'w')
    for grid in structs:
        grid = grid.lower()
        rmsds[grid] = {}
        for ligand in structs:
            ligand = ligand.lower() + '_ligand'
            print references(ligand), glides(name(ligand, grid))
            try:
                rmsds[grid][ligand] = RMSD.create(references(ligand), glides(name(ligand, grid)))
            except IOError:
                rmsds[grid][ligand] = RMSD([10.0])
            out.write(name(ligand, grid) + ':' + str(rmsds[grid][ligand]) + '\n')
elif sys.argv[1] == '-r':
    for line in open(sys.argv[2]):
        n, data = line.strip().split(':')
        ligand, grid = n.split('-to-')
        if grid not in rmsds: rmsds[grid] = {}
        rmsds[grid][ligand] = RMSD.read(data)
else:
    print 'RTFM'
    exit()


import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np


A = np.zeros((len(structs), len(structs)))
for i, grid in enumerate(structs):
    grid = grid.lower()
    for j, ligand in enumerate(structs):
        ligand = ligand.lower() + '_ligand'
        A[i, j] = rmsds[grid][ligand].get_rmsd(0)


print A

#A *= 1 / A.max()
fig, ax = plt.subplots()

heatmap = ax.pcolor(A, vmin=0, vmax=7)
plt.colorbar(heatmap)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(A.shape[1]) + 0.5, minor=False)
ax.set_yticks(np.arange(A.shape[0]) + 0.5, minor=False)
ax.invert_yaxis()
ax.xaxis.tick_top()

#lebels
column_labels = structs
row_labels = structs
ax.set_xticklabels(column_labels, minor=False, rotation = 'vertical')
ax.set_yticklabels(row_labels, minor=False)

ax.plot(ax.get_xlim(), ax.get_ylim()[::-1], linewidth = 4, c="m")
plt.savefig('top.png')
