import os
import sys

from parse_chembl import load_chembl_proc

XGLIDE_IN = '''GRIDFILE   ../../grids/{}/{}.zip
LIGANDFILE   ../../../{}/{}.mae
DOCKING_METHOD   confgen
CANONICALIZE   True
EXPANDED_SAMPLING   {}
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
WRITEREPT   True
PRECISION   {}
NENHANCED_SAMPLING   {}
'''

#ref_lig = '''USE_REF_LIGAND   True
#REF_LIGAND_FILE   ../{}
#CORE_DEFINITION   allheavy
#'''

peptide_settings = '''MAXKEEP 100000
MAXREF 1000'''

core_def = '''USE_REF_LIGAND   True
REF_LIGAND_FILE   ../../core/{}
CORE_DEFINITION   allheavy
CORE_RESTRAIN   True
CORE_POS_MAX_RMSD   2.0
CORECONS_FALLBACK   True'''

settings = { # peptides, pr, exp, enh, lig_source
    'glide1': (False, 'SP', 'True',  1, 'lig_raw'),
    'glide2': (False, 'XP', 'True',  1, 'lig_raw'),
    'glide3': (False, 'SP', 'False', 1, 'ligands/unique'),
    'glide4': (False, 'XP', 'False', 1, 'lig_raw'),
    'glide5': (False, 'SP', 'False', 1, 'lig_opt'),
    'glide6': (False, 'XP', 'False', 1, 'lig_opt'),
    'glide7': (True,  'SP', 'False', 1, 'lig_raw'),
    'glide8': (True,  'SP', 'True',  1, 'lig_raw'),
    'glide9': (True,  'SP', 'True',  4, 'lig_raw'),
    'glide10':(True,  'SP', 'False', 4, 'lig_raw'),
    'glide11':(False, 'SP', 'False', 4, 'ligands/unique'),
    'glide12':(False, 'SP', 'False', 2, 'ligands/unique')
}

queue = 'owners'

def glide_exists(ligand, grid, out_dir):
    pv_exists = os.path.exists('docking/{}/{}-to-{}/{}-to-{}_pv.maegz'.format(out_dir, ligand, grid, ligand, grid))
    rmsd_exists = os.path.exists('docking/{}/{}-to-{}/rmsd.csv'.format(out_dir, ligand, grid))
    if ligand[:6] == 'CHEMBL': return pv_exists
    return pv_exists and rmsd_exists

def glide_failed(ligand, grid, out_dir):
    log_file = 'docking/{}/{}-to-{}/{}-to-{}.log'.format(out_dir, ligand, grid, ligand, grid)
    txt1 = 'Total elapsed time'
    txt2 = 'GLIDE FATAL ERROR'
    return os.path.exists(log_file) and (txt1 in open(log_file).read() or txt2 in open(log_file).read())

def dock_pair(ligand, grid, out_dir, restrain=None):
    pv_exists = os.path.exists('docking/{}/{}-to-{}/{}-to-{}_pv.maegz'.format(out_dir, ligand, grid, ligand, grid))
    rmsd_exists = os.path.exists('docking/{}/{}-to-{}/rmsd.csv'.format(out_dir, ligand, grid))
    #if pv_exists: print ligand
    if pv_exists and (rmsd_exists or len(ligand.split('_')[0]) > 4):
        return
    if not pv_exists and glide_failed(ligand, grid, out_dir):
        return

    #if glide_exists(ligand, grid, out_dir) or glide_failed(ligand, grid, out_dir):
    #    return
    if not os.path.exists('docking/grids/{}/{}.zip'.format(grid, grid)): return
    print ligand, grid, out_dir
    os.system('mkdir -p docking/{}'.format(out_dir))
    os.chdir('docking/{}'.format(out_dir))
    time = '01:30:00'
    if '{}-to-{}'.format(ligand, grid) in os.listdir('.'):
        time = '02:30:00'
        os.system('rm -rf {}-to-{}'.format(ligand, grid))

    peptides, precision, exp, nenhanced, lig_source = settings[out_dir]

    os.system('mkdir {}-to-{}'.format(ligand, grid))
    ref_lig_file = '../../../structures/ligands/{}.mae'.format(ligand)

    with open('{}-to-{}/{}-to-{}.in'.format(ligand, grid, ligand, grid), 'w') as f:
        f.write(XGLIDE_IN.format(grid, grid, lig_source, ligand, exp, precision, nenhanced, lig_source, ligand))
        if peptides:
            f.write(peptide_settings)
        #elif ligand[:6] != 'CHEMBL': 
        #    os.chdir('../..')
        #    return # wait until the reference ligand is available

        if restrain is not None:
            f.write(core_def.format(restrain))

    with open('{}-to-{}/{}-{}.sh'.format(ligand, grid, ligand, grid), 'w') as f:
        f.write('#!/bin/bash\nmodule load schrodinger/2017-3\n')
        f.write('$SCHRODINGER/glide -WAIT {}-to-{}.in\n'.format(ligand, grid))
        if ligand[:6] != 'CHEMBL':
            f.write('$SCHRODINGER/run rmsd.py -use_neutral_scaffold -pv second -c rmsd.csv {} {}-to-{}_pv.maegz'.format(ref_lig_file, ligand, grid))

    os.chdir('{}-to-{}'.format(ligand, grid))
    os.system('sbatch -p {} -t {} -o dock.out {}-{}.sh'.format(queue, time, ligand, grid))
    os.chdir('../../..')

def dock(grids, chembl=None):
    if not os.path.exists('ligands/unique'): return
    all_ligs = sorted([l.split('.')[0] for l in os.listdir('ligands/unique')])
    to_dock = set([l for l in all_ligs if l[:6] != 'CHEMBL'])

    if chembl is not None:
        for q, c_list in chembl.items():
            to_dock.update(c_list)

    for lig in to_dock:
        if len(lig.strip()) == 0: continue
        for grid in grids:
            dock_pair(lig, grid, 'glide12')
            




