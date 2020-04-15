import os
from utils import grouper

queue = 'rondror'
group_size = 10

CMD = '{}/main.py ifp {} {} {} {}\n'

def compute_ifp(lm):
    cwd = os.getcwd()
    os.makedirs(lm.path('IFP_ROOT'), exist_ok=True)
    os.chdir(lm.path('IFP_ROOT'))
    
    ligands = lm.get_pdb() + lm.get_chembl()
    unfinished = []
    for ligand in lm.docked(ligands):
        if not os.path.exists(lm.path('IFP', {'ligand': ligand})):
            unfinished.append(ligand)
    _compute_ifp(lm, unfinished)

    os.chdir(cwd)

def _compute_ifp(lm, fp_list):
    if len(fp_list) > 0:
        print(len(fp_list), 'fp left')
    for i, ligands in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')
            for ligand in ligands:
                input_file = lm.path('DOCK_PV', {'ligand': ligand})
                output_file = lm.path('IFP', {'ligand': ligand})
                f.write(CMD.format(lm.path('CODE'), lm.params['ifp_version'],
                                   input_file, output_file,
                                   lm.params['max_poses']))
            f.write('wait\n')
        os.system('sbatch -t 02:00:00 -p {} {}fp.sh'.format(queue, i))
