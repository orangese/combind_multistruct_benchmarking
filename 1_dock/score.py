import os
import sys
from grouper import grouper
from shared_paths import shared_paths

group_size=5

output_dir = 'scores/mcss_only'
cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:}'.format(shared_paths['code'], '{}')

settings = {
    'k_list' : ['mcss'],
    'num_stats_ligs' : 10,
    'normalize' : True,
    'num_pred_chembl' : 20,
    'num_poses' : 100,
    't' : .05,
    'chembl_file': 'best_mcss.txt',
    'score_mode': 'ALL'
    #'use_chembl':False
}

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score(lm, helpers):
    all_p = [d for d in sorted(os.listdir(shared_paths['data'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot and p != 'D2R']
    settings['shared_paths'] = shared_paths
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    #os.system('rm -f *')
    write_settings_file('settings.py', settings)

    unfinished = sorted([l for l in helpers[settings['chembl_file']] 
        if not os.path.exists('{}-to-{}.sc'.format(l,lm.st))])

    if len(unfinished) > 0:
        print(len(unfinished), 'scores left')

    for i,group in enumerate(grouper(group_size, unfinished)):
        with open('{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(lm.st, lm.prot, ' '.join([q for q in group if q is not None])))
        os.system('sbatch -t 1:00:00 -p owners {}.sh'.format(i))
    os.chdir('../..')
