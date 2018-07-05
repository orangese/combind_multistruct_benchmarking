import os
import sys
from grouper import grouper

group_size=1
output_dir = 'scores/scores11'
cmd = '$SCHRODINGER/run /scratch/PI/rondror/jpaggi/combind/combind/3_analyze/scores.py {} {} {}'

settings = {
    'k_list' : ['mcss','hbond','sb2','contact'],#,'pipi'
    'num_stats_ligs' : 10,
    'normalize' : True,
    'num_pred_chembl' : 10,
    'num_poses' : 100,
    't' : 0.1,
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
    all_p = [d for d in sorted(os.listdir(lm.sp['data'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot and p != 'D2R']
    settings['shared_paths'] = lm.sp
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    #os.system('rm -f *')
    write_settings_file('settings.py', settings)

    unfinished = sorted([l for l in helpers[settings['chembl_file']] 
        if not os.path.exists('{}-to-{}.sc'.format(l,lm.st))])

    if len(unfinished) > 0:
        print len(unfinished), 'scores left'

    for i,group in enumerate(grouper(group_size, unfinished)):
        with open('{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(lm.st, lm.prot, ' '.join([q for q in group if q is not None])))
        os.system('sbatch -t 1:00:00 -p rondror {}.sh'.format(i))
    os.chdir('../..')
