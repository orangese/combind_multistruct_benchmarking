##!/share/PI/rondror/software/schrodinger2016-1/run

import os
import subprocess
from multiprocessing import Pool
from tests import get_first
GLIDE = 'glide'#'/share/PI/rondror/software/schrodinger2017-1/glide'
#GLIDE = '/share/software/modules/chemistry/schrodinger2017-2/glide'
XGLIDE_IN = '''GRIDFILE   ../../grids/{}/{}.zip
LIGANDFILE   ../../unique_ligands/{}.mae
USE_REF_LIGAND   True
REF_LIGAND_FILE   ../../unique_ligands/{}.mae
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

def glide_exists(ligand, grid):
    return os.path.exists('{}-to-{}/{}-to-{}_pv.maegz'.format(ligand, grid, ligand, grid))

def glide_failed(ligand, grid):
    log_file = '{}-to-{}/{}-to-{}.log'.format(ligand, grid, ligand, grid)
    return os.path.exists(log_file) and 'Total elapsed time' in open(log_file).read()

def dock(pair):
    ligand, grid = pair
    #print 'docking {} to {}'.format(ligand, grid)
    if '{}-to-{}'.format(ligand, grid) in os.listdir('.'):
        os.system('rm -rf {}-to-{}'.format(ligand, grid))

    os.system('mkdir {}-to-{}'.format(ligand, grid))
    with open('{}-to-{}/{}-to-{}.in'.format(ligand, grid, ligand, grid), 'w+') as f:
        f.write(XGLIDE_IN.format(grid, grid, ligand, ligand))

    os.chdir('{}-to-{}'.format(ligand, grid))

    subprocess.call([GLIDE, "{}-to-{}.in".format(ligand, grid), "-WAIT"])
    os.chdir('..')
    #print '{} to {} returned'.format(ligand, grid)
    return (ligand, grid, glide_exists(ligand, grid) or glide_failed(ligand, grid))

def dock_dataset(print_only=False):
    os.system('mkdir -p xglide')
    assert len(os.listdir('grids')) == 1
    grid = os.listdir('grids')[0]
    os.chdir('xglide')

    to_dock = []

    all_ligands = [f.split('.')[0] for f in os.listdir('../unique_ligands')]
    for ligand in all_ligands:
        if not glide_exists(ligand, grid) and not glide_failed(ligand, grid):
            to_dock.append((ligand, grid))
    
    if len(to_dock) > 0: print grid, len(to_dock)
    if print_only: return
    #print to_dock
    #return
    num_licenses = 10
    pool = Pool(num_licenses)

    i = 0
    while len(to_dock) > 0 and i < 5:
        i += 1
        print 'iteration {}, {} jobs left to go'.format(i, len(to_dock))

        currently_docking = to_dock
        to_dock = []

        done = 0
        not_done = 0
        for l, g, success in pool.imap_unordered(dock, currently_docking):
            if not success:
                to_dock.append((l, g))
                not_done += 1
            else: done += 1
        print '{} done, {} not done'.format(done, not_done)
        if len(to_dock) == 0:
            break

    else:
        print len(to_dock), 'failed pairs.'
        print to_dock
