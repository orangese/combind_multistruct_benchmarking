#!/share/PI/rondror/software/schrodinger2016-1/run

import os
import slurm
from multiprocessing import Pool

GLIDE = '/share/PI/rondror/software/schrodinger2017-1/glide'

XGLIDE_IN = '''GRIDFILE   ../../grids/{}/{}.zip
LIGANDFILE   ../../ligands/{}.mae
USE_REF_LIGAND   True
REF_LIGAND_FILE   ../../ligands/{}.mae
CORE_DEFINITION   allheavy
DOCKING_METHOD   confgen
EXPANDED_SAMPLING   True
NENHANCED_SAMPLING   4
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
NMAXRMSSYM   300
MAXKEEP   100000
MAXREF   1000
PRECISION   SP
WRITEREPT   True
'''

def glideExists(ligand, grid):
    return os.path.exists('{}-to-{}/{}-to-{}_pv.maegz'.format(ligand, grid, ligand, grid)) 

def glideFailed(ligand, grid):
    logFile = '{}-to-{}/{}-to-{}.log'.format(ligand, grid, ligand, grid)
    return os.path.exists(logFile) and 'Total elapsed time' in open(logFile).read()

def dock(pair):
    ligand, grid = pair

    if '{}-to-{}'.format(ligand, grid) in os.listdir('.'):
        os.system('rm -rf {}-to-{}'.format(ligand, grid))

    os.system('mkdir {}-to-{}'.format(ligand, grid))
    
    os.chdir('{}-to-{}'.format(ligand, grid))
    with open('{}-to-{}.in'.format(ligand, grid), 'w+') as f:
        f.write(XGLIDE_IN.format(grid, grid, ligand, ligand))

    slurm.salloc('{} {}-to-{}.in -WAIT'.format(GLIDE, ligand, grid), 1, "02:30:00")
    os.chdir('..')

    return (ligand, grid, glideExists(ligand, grid) or glideFailed(ligand, grid))

def dockDataset():
    os.system('mkdir -p xglide')
    os.chdir('xglide')
        
    toDock = []

    for ligand in os.listdir('../ligands'):
        ligand = ligand.split('.')[0]
        for grid in os.listdir('../grids'):
            if not glideExists(ligand, grid) and not glideFailed(ligand, grid):
                toDock.append((ligand, grid))

    num_licenses = 5
    pool = Pool(num_licenses)

    i = 0    
    while len(toDock) > 0:
        i += 1
        print 'iteration {}, {} jobs left to go'.format(i, len(toDock))

        currentlyDocking = toDock
        toDock = []
        
        done = 0 
        not_done = 0
        for l, g, success in pool.imap_unordered(dock, currentlyDocking):
            if not success:
                toDock.append((l, g))
                not_done += 1
            else: done += 1
        print '{} done, {} not done'.format(done, not_done)
        if len(toDock) == 0:
            print 'all done!'
            break

    else:
        print len(toDock), 'failed pairs.'
        print toDock
