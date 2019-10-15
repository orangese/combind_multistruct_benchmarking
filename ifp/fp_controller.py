import os
import sys
from utils import grouper
from settings import paths

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/main.py ifp -version {} -mode pv -input_file {} -output_file {} -poses {}\n'
st_cmd = '$SCHRODINGER/run {}/main.py ifp -version {} -mode st -output_file {}\n'

def get_fp(lm, fp_list):
    if len(fp_list) > 0:
        print(len(fp_list), 'fp left')
    for i, names in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for name in names:
                if name is None: continue
                input_file = lm.path('PV', {'ligand': name})
                output_file = lm.path('IFP', {'ligand': name})
                f.write(pv_cmd.format(lm.path('CODE'), input_file, output_file,
                                      lm.params['max_poses']))
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
            f.write(st_cmd.format(paths['CODE'], output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def compute_fp(lm, pdb=False):
    os.system('mkdir -p ifp')
    os.system('mkdir -p ifp/{}'.format(lm.path['IFP_ROOT']))
    ligands = lm.pdb + lm.chembl()
    unfinished = []
    for ligand in lm.docked(ligands):
        if pdb and 'CHEMBL' in lig: continue
        if not os.path.exists(lm.path('IFP', {'ligand': name})):
            unfinished.append(ligand)
    
    os.chdir(lm.path('IFP_ROOT')) 
    structure_fp(lm)
    get_fp(lm, unfinished, raw)
    os.chdir(lm.path('ROOT'))

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
