import os
import sys
from grouper import grouper
from shared_paths import shared_paths

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/main.py ifp -mode pv -input_file {} -output_file {} -raw {} -write_root {} -read_root {}\n'
st_cmd = '$SCHRODINGER/run {}/main.py ifp -mode st -output_file {} -write_root {} -read_root {}\n'

def get_fp(lm, fp_list, raw):
    if len(fp_list) > 0:
        print(len(fp_list), 'fp left')
    for i, names in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for name in names:
                if name is None: continue
                input_file = '../../docking/{}/{}/{}_pv.maegz'.format(shared_paths['docking'],
                                                                      name, name)
                output_file = '{}-{}.fp'.format(name, shared_paths['docking'])
                f.write(pv_cmd.format(shared_paths['code'], input_file, output_file, raw, lm.write_root, lm.read_root)) 
            f.write('wait\n')
        os.system('sbatch --cpus-per-task=1 --time=02:00:00 -p {} {}fp.sh'.format(queue, i))

def structure_fp(lm):
    for lig in sorted(os.listdir('{}/structures/ligands'.format(lm.read_root))):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print ('structure fp', pdb)
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(st_cmd.format(shared_paths['code'], output_file, lm.write_root, lm.read_root)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def compute_fp(lm, raw = False):
    os.system('mkdir -p ifp')
    os.system('mkdir -p ifp/{}'.format(shared_paths['ifp']['version']))
    ligands = lm.pdb
    if not raw: ligands += lm.chembl()
    unfinished = []
    for lig in lm.docked(ligands):
        name = '{}-to-{}'.format(lig, lm.st, shared_paths['docking'])
        output_file = 'ifp/{}/{}-{}.fp'.format(shared_paths['ifp']['version'],
                                               name, shared_paths['docking'])
        if not os.path.exists(output_file):
            unfinished.append(name)
    
    os.chdir('ifp/{}'.format(shared_paths['ifp']['version'])) 
    if not raw: structure_fp(lm)
    get_fp(lm, unfinished, raw)
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
