#from rmsd import RMSD
import sys
import os

#os.chdir('/share/PI/rondror/docking_code/3_score/')
from rmsd import RMSD

grids_dir = sys.argv[1]# '/scratch/PI/rondror/docking_data/B2AR/grids'
OUTPUT = 'rmsd.csv'
os.chdir(grids_dir)
os.chdir('..')
structs = [d for d in os.listdir(grids_dir) if os.path.isdir(os.path.join(grids_dir, d))]
print(structs)

def references(ligand):
    return "ligands/{}.mae".format(ligand)

def name(ligand, grid):
    return "{}-to-{}".format(ligand, grid)

def glides(name):
    return "glide/{}/{}_pv.maegz".format(name, name)
print(os.system('ls'))    
rmsds = {}
if True:#sys.argv[1] == '-o':
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
elif sys.argv[1] == '-r':
    for line in open(sys.argv[2]):
        n, data = line.strip().split(':')
        ligand, grid = n.split('-to-')
        if grid not in rmsds: rmsds[grid] = {}
        rmsds[grid][ligand] = RMSD.read(data)
else:
    print 'RTFM'
    exit()
