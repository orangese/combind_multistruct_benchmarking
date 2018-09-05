import os
import sys
from shared_paths import shared_paths
sys.path.append('../3_analyze')
from containers import LigandManager


output_dir = 'scores/pdb'
cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:}'.format(shared_paths['code'], '{}')

settings = {
    'k_list' : ['mcss', 'sb2', 'contact', 'hbond', 'pipi'],
    'num_stats_ligs' : 10,
    'normalize' : True,
    'num_poses' : 100,
    't' : .1,
    'score_mode': 'ALL',
    'chembl': False
}

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score_pdb(lm):
    all_p = [d for d in sorted(os.listdir(lm.sp['data'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot and p != 'D2R']
    settings['shared_paths'] = lm.sp
    write_settings_file('settings.py', settings)

    with open('pdb.sh','w') as f:
        f.write('#!/bin/bash\n')
        f.write(cmd.format(lm.st, lm.prot, ' '.join(lm.docked(lm.pdb)[:10]))+'\n')
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
