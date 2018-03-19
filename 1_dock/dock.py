import os
import sys

import core_proc as core
from parse_chembl import load_chembl_proc, load_drugs

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
    'glide3': (False, 'SP', 'False', 1, 'lig_raw'),
    'glide4': (False, 'XP', 'False', 1, 'lig_raw'),
    'glide5': (False, 'SP', 'False', 1, 'lig_opt'),
    'glide6': (False, 'XP', 'False', 1, 'lig_opt'),
    'glide7': (True,  'SP', 'False', 1, 'lig_raw'),
    'glide8': (True,  'SP', 'True',  1, 'lig_raw'),
    'glide9': (True,  'SP', 'True',  4, 'lig_raw'),
    'glide10':(True,  'SP', 'False', 4, 'lig_raw'),
    'glide11':(False, 'SP', 'False', 4, 'lig_raw'),
    'glide12':(False, 'SP', 'False', 2, 'ligands/unique')
}

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

def dock(ligand, grid, out_dir=None, restrain=None):
    if out_dir is None: out_dir = 'glide12'
    pv_exists = os.path.exists('docking/{}/{}-to-{}/{}-to-{}_pv.maegz'.format(out_dir, ligand, grid, ligand, grid))
    rmsd_exists = os.path.exists('docking/{}/{}-to-{}/rmsd.csv'.format(out_dir, ligand, grid))
    #if pv_exists: print ligand
    if pv_exists and (rmsd_exists or ligand[:6] == 'CHEMBL'):
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

    peptides, precision, exp, nenhanced, lig_source = settings['glide12']#out_dir]

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
    os.system('sbatch -p rondror -t {} -o dock.out {}-{}.sh'.format(time, ligand, grid))
    os.chdir('../../..')

#from parse_chembl import load_chembl_proc

def dock_dataset(grids, num_chembl=0):

    chembl_info = load_chembl_proc()
    drugs = load_drugs()
    all_ligs = sorted([l.split('_')[0] for l in os.listdir('ligands/unique')])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL'] + [l for l in all_ligs if l in drugs]
    chembl_ligs = [l for l in all_ligs if l in chembl_info and chembl_info[l].valid_stereo]
    chembl_ligs.sort(key=lambda x: chembl_info[x].ki)

    for lig in pdb_ligs + chembl_ligs[:num_chembl]:
        for grid in grids:#grids:#sorted(os.listdir('docking/grids')):
            dock('{}_lig'.format(lig), grid, 'glide12')
            
def core_dock(grids):
    matches = core.load_matches()
    for ss in matches:
        #if ss != '6C38ss.mae': continue
        for lig in matches[ss]:
            if not matches[ss][lig]: continue
            for grid in grids:
                dock(lig, grid, 'glide12-{}'.format(ss.split('.')[0]), ss)
                #break




