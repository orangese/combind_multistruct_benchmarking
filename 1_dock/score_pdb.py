import os
import sys
from shared_paths import shared_paths
from containers import LigandManager

max_ligands = 20

output_dir = 'scores/pdb_seperated_hbond_pipi/'
cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:}'.format(shared_paths['code'], '{}')

settings = {
    'k_list' : ['mcss', 'sb2', 'pipi', 'contact', 'hbond_donor', 'hbond_acceptor'],
    'num_stats_ligs' : 10,
    'num_poses' : 100,
    'chembl': False
}

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score_pdb(lm):
    ligands = lm.docked(lm.pdb)[:max_ligands+1]
    self_docked = lm.st+'_lig'
    if self_docked in ligands:
        ligands.remove(self_docked)
    else:
        ligands.pop(-1)
    if len(ligands) == 1: return
    all_p = [d for d in sorted(os.listdir(lm.sp['data'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot and p != 'D2R']
    settings['shared_paths'] = lm.sp
    settings['t'] = 1 / float(len(ligands)-1)
    write_settings_file('settings.py', settings)

    with open('pdb.sh','w') as f:
        f.write('#!/bin/bash\n')
        f.write(cmd.format(lm.st, lm.prot, ' '.join(ligands))+'\n')
    os.system('sbatch -t 1:00:00 -p owners pdb.sh')

datasets = sys.argv[1:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir(shared_paths['data'])) if d[0] != '.' and d[-3:] != 'old']

os.chdir(shared_paths['data'])
for i, d in enumerate(datasets):
    print(d, i)
    os.chdir(d)
    lm = LigandManager(shared_paths, d)
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    score_pdb(lm)
    os.chdir('../../..')
