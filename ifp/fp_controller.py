import os
import sys
from grouper import grouper
from shared_paths import shared_paths

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/main.py ifp {} -mode pv -input_file {} -output_file {} -raw {}\n'
st_cmd = '$SCHRODINGER/run {}/main.py ifp {} -mode st -output_file {}\n'

def get_fp(protein, lm, fp_list, raw):
    ''' Create fingerprinting batch scripts for all docked ligands.
    Note: fp_list == unfinished
    '''
    if len(fp_list) > 0:
        print('Fingerprinting {} ligands'.format(len(fp_list)))
    for i, names in enumerate(grouper(group_size, fp_list)):
        with open('fp-group{}.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for name in names:
                if name is None: continue
                input_file = '../../docking/{}/{}/{}_pv.maegz'.format(shared_paths['docking'],
                                                                      name, name)
                output_file = '{}-{}.fp'.format(name, shared_paths['docking'])
                f.write(pv_cmd.format(shared_paths['code'], protein, input_file, output_file, raw)) 
            f.write('wait\n')
        os.system('sbatch --cpus-per-task=1 --time=02:00:00 -p {} fp-group{}.sh'.format(queue, i))

def structure_fp(protein, lm):
    ''' For each PDB ligand without a structure fingerprint computed already, compute it now.
    '''
    for lig in sorted(os.listdir('{}{}/structures/ligands'.format(shared_paths['read_data'],protein))):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print ('Create structure fp for ligand from ', pdb)
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(st_cmd.format(shared_paths['code'], protein, output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def compute_fp(protein, lm, raw = False):
    ''' Must be called from write_data/

    Inputs:
    * raw (bool) (optionali, default=False): if raw is set to True, this script will compute fingerprints
    for chembl and PDB ligands. Else, the script will only compute fingerprints for PDB ligands
    '''
    os.system('mkdir -p ifp')
    os.system('mkdir -p ifp/{}'.format(shared_paths['ifp']['version']))

    # Choose ligands to use depending on raw
    ligands = lm.pdb
    if not raw: ligands += lm.chembl()

    # For each ligand, check for its fingerprint output file. This file should take the form
    # "ligid-to-grid-dockingversion.fp". If output has not been generated, append the ligand to unfinished
    unfinished = []
    for lig in lm.docked(ligands):
        name = '{}-to-{}'.format(lig, lm.st)
        output_file = 'ifp/{}/{}-{}.fp'.format(shared_paths['ifp']['version'], name, shared_paths['docking'])
        if not os.path.exists(output_file):
            unfinished.append(name)
    
    os.chdir('ifp/{}'.format(shared_paths['ifp']['version'])) 

    if not raw: structure_fp(protein, lm)
    get_fp(protein, lm, unfinished, raw)

    os.chdir('../..')

def parse_fp_file(fp_file):
    ifps = {}
    try:
        with open(fp_file) as f:
            pose_num = 0
            for line in f:
                if line.strip() == '': continue
                if line[:4] == 'Pose':
                    pose_num = int(line.strip().split(' ')[1])
                    ifps[pose_num] = {}
                    continue
                sc_key, sc = line.strip().split('=')
                i,r,ss = sc_key.split('-')
                i = int(i)
                sc = float(sc)
                prev_sc = ifps[(i, r)] if (i,r) in ifps[pose_num] else 0
                ifps[pose_num][(i,r)] = max(prev_sc, sc)

    except Exception as e:
        print(e)
        print(fp_file, 'fp not found')
    if len(ifps) == 0:
        print('check', fp_file)
        return {}
    return ifps
