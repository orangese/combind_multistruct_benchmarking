import os
import sys

from shared_paths import shared_paths
from grouper import grouper
from pick_helpers import load_helpers
sys.path.append('../3_analyze')
from containers import LigandManager

group_size=30
cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:}'.format(shared_paths['code'], '{}')

def write_settings_file(out_path, settings):
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score(lm, helpers, settings, subdir):
    all_p = [d for d in sorted(os.listdir(lm.sp['data'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot and p != 'D2R']
    settings['shared_paths'] = lm.sp
    os.system('mkdir -p {}'.format(subdir))
    os.chdir(subdir)
    write_settings_file('settings.py', settings)
    print(os.getcwd())
    unfinished = sorted([l for l in helpers[settings['chembl_file']]
                         if not os.path.exists('{}-to-{}.sc'.format(l,lm.st))])

    if len(unfinished) > 0:
        print(len(unfinished), 'scores left')

    for i,group in enumerate(grouper(group_size, unfinished)):
        with open('{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(lm.st, lm.prot, ' '.join([q for q in group if q is not None])))
        os.system('sbatch -t 1:00:00 -p owners {}.sh'.format(i))
    os.chdir('..')

#############################################################################

output_dir, chembl_file, features, temps, n_ligs = sys.argv[1:6]
features = features.split(',')
temps = [float(t) for t in  temps.split(',')]
n_ligs = [int(n) for n in n_ligs.split(',')]

datasets = sys.argv[6:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir(shared_paths['data'])) if d[0] != '.' and d[-3:] != 'old']

os.chdir(shared_paths['data'])
for i, d in enumerate(datasets):
    print(d, i)
    os.chdir(d)
    lm = LigandManager(shared_paths, d)
    helpers = load_helpers()
    os.chdir('scores')
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    for t in temps:
        print(t)
        for n in n_ligs:
            subdir = "t={},n={}".format(t, n)
            settings = {
                'num_pred_chembl' : n,
                't' : t,
                'k_list' : features,
                'chembl_file': chembl_file,
                'num_stats_ligs' : 10,
                'normalize' : True,
                'num_poses' : 100,
                'score_mode': 'ALL'
            }
            score(lm, helpers, settings, subdir)
    os.chdir('../../..')
