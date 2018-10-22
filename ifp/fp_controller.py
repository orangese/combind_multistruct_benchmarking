import os
import sys
from grouper import grouper
from shared_paths import shared_paths

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/ifp/fp.py -mode pv -input_file {} -output_file {}\n'
st_cmd = '$SCHRODINGER/run {}/ifp/fp.py -mode st -output_file {}\n'

def get_fp(lm, fp_list):
    if len(fp_list) > 0:
        print(len(fp_list), 'fp left')
    for i, pairs in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for p in pairs:
                if p is None: continue
                input_file = '../../docking/{}/{}/{}_pv.maegz'.format(shared_paths['docking'], p, p)
                output_file = '{}.fp'.format(p)
                f.write(pv_cmd.format(shared_paths['code'], input_file, output_file)) 
            f.write('wait\n')
        os.system('sbatch --cpus-per-task=1 --time=02:00:00 -p {} {}fp.sh'.format(queue, i))

def structure_fp(lm):
    for lig in sorted(os.listdir('../../structures/ligands')):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print ('structure fp', pdb)
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(st_cmd.format(shared_paths['code'], output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def compute_fp(lm):
    os.system('mkdir -p ifp')
    os.system('mkdir -p ifp/{}'.format(shared_paths['ifp']))

    unfinished = []    
    for lig in lm.docked(lm.pdb+lm.chembl()):
        output_file = 'ifp/{}/{}-to-{}.fp'.format(shared_paths['ifp'], lig, lm.st)
        if os.path.exists(output_file): continue
        unfinished.append('{}-to-{}'.format(lig, lm.st))
       
    os.chdir('ifp/{}'.format(shared_paths['ifp'])) 
    structure_fp(lm) # Should move this so it is called by itself.
    get_fp(lm, unfinished)
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
