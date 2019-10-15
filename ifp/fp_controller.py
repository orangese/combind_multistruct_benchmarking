import os
import sys
from utils import grouper
from settings import paths

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/main.py ifp -version {} -mode pv -input_file {} -output_file {} -poses {}\n'
st_cmd = '$SCHRODINGER/run {}/main.py ifp -version {} -mode st -output_file {}\n'

def compute_fp(lm, pdb=False):
    root = lm.path('IFP_ROOT')
    os.system('mkdir -p {}'.format('/'.join(str(root).split('/')[:-1])))
    os.system('mkdir -p {}'.format(root))

    cwd = os.getcwd()
    os.chdir(lm.path('IFP_ROOT'))

    _structure_fp(lm)
    
    ligands = lm.pdb + lm.chembl()
    unfinished = []
    for ligand in lm.docked(ligands):
        if pdb and 'CHEMBL' in lig: continue
        if not os.path.exists(lm.path('IFP', {'ligand': ligand})):
            unfinished.append(ligand)
    _compute_fp(lm, unfinished)

    os.chdir(cwd)

def _structure_fp(lm):
    for ligand in sorted(os.listdir('../../structures/ligands')):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print ('structure fp', pdb)
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(st_cmd.format(paths['CODE'], output_file))
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def _compute_fp(lm, fp_list):
    if len(fp_list) > 0:
        print(len(fp_list), 'fp left')
    for i, ligands in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for ligand in ligands:
                input_file = lm.path('DOCK_PV', {'ligand': ligand})
                output_file = lm.path('IFP', {'ligand': ligand})
                f.write(pv_cmd.format(lm.path('CODE'), lm.params['ifp_version'],
                                      input_file, output_file,
                                      lm.params['max_poses']))
            f.write('wait\n')
        os.system('sbatch -t 02:00:00 -p {} {}fp.sh'.format(queue, i))

################################################################################
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
