import os
import sys

from shared_paths import shared_paths, proteins
from grouper import grouper
from pick_helpers import load_helpers
from containers import Protein

group_size=30
cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:}'.format(shared_paths['code'], '{}')

def write_settings_file(out_path, settings):
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score(lm, helpers, settings, subdir):
    settings['stats_prots'] = [p for p in proteins if p != lm.protein]
    settings['shared_paths'] = shared_paths
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
            f.write(cmd.format(lm.st, lm.protein, ' '.join([q for q in group if q is not None])))
        os.system('sbatch -t 1:00:00 -p owners {}.sh'.format(i))
    os.chdir('..')

#############################################################################

output_dir, chembl_file, features, n_ligs = sys.argv[1:5]
features = features.split(',')
n_ligs = [int(n) for n in n_ligs.split(',')]

datasets = sys.argv[5:]
if datasets == []:
    datasets = proteins

os.chdir(shared_paths['data'])
for i, d in enumerate(datasets):
    print(d, i)
    os.chdir(d)
    prot = Protein(d)
    lm = prot.lm

    helpers = load_helpers()
    self_docked = lm.st+'_lig'
    if self_docked in helpers:
        del helpers[self_docked]
    os.chdir('scores')
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    for n in n_ligs:
        subdir = "n={}".format(n)
        settings = {
            'num_pred_chembl' : n,
            't' : 0.8 / float(n),
            'k_list' : features,
            'chembl_file': chembl_file,
            'num_stats_ligs' : 20,
            'num_poses' : 100,
            'chembl': True
        }
        score(lm, helpers, settings, subdir)
    os.chdir('../../..')
