import os
import sys
from grouper import grouper
from shared_paths import shared_paths
from schrodinger.structure import StructureReader

XGLIDE_IN = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
EXPANDED_SAMPLING   False
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
WRITEREPT   True
PRECISION   SP
NENHANCED_SAMPLING   2
'''

GLIDE_ES4 = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
EXPANDED_SAMPLING   False
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
WRITEREPT   True
PRECISION   SP
NENHANCED_SAMPLING   4
'''

GLIDE_ES1 = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
EXPANDED_SAMPLING   False
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
WRITEREPT   True
PRECISION   SP
NENHANCED_SAMPLING   1
'''

GLIDEXP = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
WRITEREPT   True
PRECISION   XP
'''

EXPANDED = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
EXPANDED_SAMPLING   True
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
NMAXRMSSYM   300
MAXKEEP   100000
MAXREF   1000
WRITEREPT   True
PRECISION   SP
NENHANCED_SAMPLING   4
'''

INPLACE = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   inplace
WRITEREPT   True
PRECISION   SP
'''

REFINE = '''GRIDFILE   {}{}/{}.zip
LIGANDFILE   ../../../ligands/prepared_ligands/{}/{}.mae
DOCKING_METHOD   mininplace
WRITEREPT   True
PRECISION   SP
'''

modes = {'confgen':{
                       'name':     shared_paths['docking'],
                       'template': XGLIDE_IN},
        'confgen_es1':{
                       'name':     'confgen_es1',
                       'template': GLIDE_ES1},
        'confgen_es4': {
                       'name':     'confgen_es4',
                       'template': GLIDE_ES4},
        'inplace':{
                       'name':     'inplace',
                       'template': INPLACE},
        'mininplace':{
                       'name':    'mininplace',
                       'template': REFINE},
        'expanded':{
                       'name':     'expanded',
                       'template': EXPANDED},
        'XP':{
                       'name':     'XP',
                       'template': GLIDEXP},
        }

queue = 'owners'
group_size = 5

dock_cmd = '$SCHRODINGER/glide -WAIT {}-to-{}.in\n' 
rmsd_cmd = '$SCHRODINGER/run rmsd.py -use_neutral_scaffold -pv second -c rmsd.csv ../../../structures/ligands/{}.mae {}-to-{}_pv.maegz\n'

def get_state(ligand, grid, grid_root):
    ''' Detect the stage of processing of this ligand/grid combination. The state is returned as an int
    between 0 and 2, inclusive. The states are as follows:
    * 0: do nothing
    * 1: compute rmsd
    * 2: dock
    '''
    pv = '{}-to-{}/{}-to-{}_pv.maegz'.format(ligand, grid, ligand, grid)
    rmsd = '{}-to-{}/rmsd.csv'.format(ligand, grid)
    log = '{}-to-{}/{}-to-{}.log'.format(ligand, grid, ligand, grid)
    ref_lig = '../../structures/ligands/{}.mae'.format(ligand)
    inp_lig = '../../ligands/prepared_ligands/{}/{}.mae'.format(ligand, ligand)
    inp_grid = '{}{}/{}.zip'.format(grid_root, grid, grid)

    if os.path.exists(rmsd): 
        return 0
    if os.path.exists(pv):
        if os.path.exists(ref_lig): return 1
        else: return 0
    if os.path.exists(log):
        logtxt = open(log).read()
        if 'Total elapsed time' in logtxt or 'GLIDE FATAL ERROR' in logtxt:
            return 0
    if not os.path.exists(inp_lig) or not os.path.exists(inp_grid):
        return 0

    os.system('rm -rf {}-to-{}'.format(ligand, grid))
    return 2

def write_inp_files(grid_root, all_pairs, mode):

    TEMPLATE = modes[mode]['template']

    os.system('rm *.sh')
    for ligand, grid in all_pairs:
        os.system('mkdir {}-to-{}'.format(ligand, grid))
        with open('{}-to-{}/{}-to-{}.in'.format(ligand, grid, ligand, grid), 'w') as f:
            f.write(TEMPLATE.format(grid_root, grid, grid, ligand, ligand))

def proc_all(all_pairs, dock=False, rmsd=False):
    
    for i,group in enumerate(grouper(group_size, all_pairs)):
        with open('dock{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            for pair in group:
                if pair is None: continue
                ligand, grid = pair
                f.write('cd {}-to-{}\n'.format(ligand, grid))
                if dock:
                    f.write(dock_cmd.format(ligand, grid))
                if rmsd:
                    f.write(rmsd_cmd.format(ligand, ligand, grid))
                f.write('cd ..\n')

        os.system('sbatch -p {} -t 1:00:00 -o dock.out dock{}.sh'.format(queue, i))

def dock(lm, protein, dude=None, maxnum=30, mode = 'confgen'):
    ''' dock() modified to only dock dude ligands

    Inputs:
    * dude (dict): maps dude_helper_filename (str) -> pdb_ligand_helpers (dict)
        * pdb_ligand_helpers (dict): maps ID of a PDB ligand -> list of IDs of all dude helpers
    '''

    read_root = '{}{}/'.format(shared_paths['read_data'], protein)
    grid_root = '{}/docking/grids/'.format(read_root)
    write_root = '{}{}/'.format(shared_paths['write_data'], protein)

    if lm.st is None: return
    docking = modes[mode]['name']
    os.system('mkdir -p docking/{}'.format(docking))
    os.chdir('docking/{}'.format(docking))

    # If not called from prepare_all, use the first maxnum pdb ligands
    if dude is None:
        ligs = lm.pdb[:maxnum]
        grids = [lm.st]
    else:
        ligs = set([])
        grids = [lm.st]
        for f, f_data in dude.items():
            for q, c_list in f_data.items():
                ligs.add(q)
                ligs.update(c_list)

    # print("Docking {} ligands on {} grids".format(len(ligs), len(grids)))
    to_dock = []
    to_rmsd = []
    for lig in ligs:
        for grid in grids:
            s = get_state(lig, grid, grid_root)
            if s == 2: to_dock.append((lig, grid))
            if s == 1: to_rmsd.append((lig, grid))
            
    if len(to_dock) > 0:
        print('Docking {} ligands'.format(len(to_dock)))
        write_inp_files(grid_root, to_dock, mode)
        proc_all(to_dock, dock=True)
    if len(to_rmsd) > 0:
        print('Computing {} rmsds'.format(len(to_rmsd)))
        proc_all(to_rmsd, rmsd=True)

    # os.chdir('../..')
    os.chdir(write_root)

def verify_dock(lm):
    for ligand in lm.docked(lm.all_ligs):
        prepared = "{0:}/ligands/prepared_ligands/{1:}/{1:}.mae".format(lm.root, ligand)
        docked = "{0:}/docking/{1:}/{2:}-to-{3:}/{2:}-to-{3:}_pv.maegz".format(lm.root,
                                                                               shared_paths['docking'],
                                                                               ligand,
                                                                               lm.st)
        
        if not os.path.exists(prepared) and os.path.exists(docked):
            print('Ligand or docking files do not exist for:')
            print(prepared)
            print(docked)
            continue

        prepared = next(StructureReader(prepared))
        docked = list(StructureReader(docked))[1]
        if not prepared.isEquivalent(docked):
            print('{} is not the same in glide output and prepared ligands.'.format(ligand))
