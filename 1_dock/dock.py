import os
import sys

import itertools

group_size=10
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

XGLIDE_IN = '''GRIDFILE   ../../grids/{}/{}.zip
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

queue = 'rondror'
out_dir = 'docking/glide12'

dock_cmd = '$SCHRODINGER/glide -WAIT {}-to-{}.in\n' 
rmsd_cmd = '$SCHRODINGER/run rmsd.py -use_neutral_scaffold -pv second -c rmsd.csv ../../../structures/ligands/{}.mae {}-to-{}_pv.maegz\n'

def get_state(ligand, grid):
    pv = '{}-to-{}/{}-to-{}_pv.maegz'.format(ligand, grid, ligand, grid)
    rmsd = '{}-to-{}/rmsd.csv'.format(ligand, grid)
    log = '{}-to-{}/{}-to-{}.log'.format(ligand, grid, ligand, grid)
    ref_lig = '../../structures/ligands/{}.mae'.format(ligand)
    inp_lig = '../../ligands/prepared_ligands/{}/{}.mae'.format(ligand, ligand)
    inp_grid = '../grids/{}/{}.zip'.format(grid, grid)

    # 0: do nothing
    # 1: compute rmsd
    # 2: dock

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

def write_inp_files(all_pairs):
    os.system('rm -f *.sh')
    for ligand, grid in all_pairs:
        os.system('mkdir {}-to-{}'.format(ligand, grid))
        with open('{}-to-{}/{}-to-{}.in'.format(ligand, grid, ligand, grid), 'w') as f:
            f.write(XGLIDE_IN.format(grid, grid, ligand, ligand))

def proc_all(all_pairs, dock=False, rmsd=False):
    
    for i,group in enumerate(grouper(group_size, all_pairs)):
        with open('dock{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
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

def dock(lm, chembl=None, maxnum=20):

    os.system('mkdir -p {}'.format(out_dir))
    os.chdir(out_dir)

    if chembl is None:
        ligs = lm.pdb[:maxnum]
        grids = lm.grids[:maxnum]
    else:
        ligs = set([])
        grids = [lm.st]
        for f, f_data in chembl.items():
            for q, c_list in f_data.items():
                ligs.add(q)
                ligs.update(c_list)

    to_dock = []
    to_rmsd = []
    for lig in ligs:
        for grid in grids:
            s = get_state(lig, grid)
            if s == 2: to_dock.append((lig, grid))
            if s == 1: to_rmsd.append((lig, grid))
            
    if len(to_dock) > 0:
        print 'docking {} pairs'.format(len(to_dock))
        write_inp_files(to_dock)
        proc_all(to_dock, dock=True)
    if len(to_rmsd) > 0:
        print 'computing {} rmsds'.format(len(to_rmsd))
        proc_all(to_rmsd, rmsd=True)

    os.chdir('../..')

