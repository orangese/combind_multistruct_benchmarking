import os
import sys

sys.path.append('../3_analyze')
from containers import Dataset

output_dir = 'scores/test'
script = '/scratch/PI/rondror/jbelk/method/combind/4_score/score_query.py'

settings = {
    'data_dir' : '/scratch/PI/rondror/jbelk/method/data',
    'glide_dir' : 'docking/glide12',
    'ifp_dir' : 'ifp/ifp2',
    'mcss_dir' : 'mcss/mcss7',
    'features' : {'mcss':[],'hbond':[2,3]},#,'pipi':[6],'contact':[11]},
    #'stats_prots' : ['5HT2B','SIGMA1'],
    #'struct_dict' : struct_dict,
    'num_stats_ligs' : 10,
    'num_stats_poses' : 25,
    'smooth' : 0.02,
    'normalize' : True,
    'num_pred_chembl' : 12,
    'num_poses' : 100,
    't' : 10,
    'mcss_sort': False, # set to True when using best_mcss.txt
    'chembl_file': 'best_affinity.txt'
    #'req_stereo':True,
    #'use_chembl':False
}

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def write_job_file(out_path, query, struct, prot):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('$SCHRODINGER/run {} {} {} {}\n'.format(script, query, struct, prot))

def score(prot, struct, helpers):
    settings['stats_prots'] = ['AR','SIGMA1','5HT2B','CHK1']
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir('{}'.format(output_dir))
    write_settings_file('settings.py', settings)
    for q in helpers[settings['chembl_file']]:
        write_job_file('{}-to-{}.in'.format(q, struct), q, struct, prot)
        os.system('sbatch -p rondror -t 1:00:00 -o {}-to-{}.out {}-to-{}.in'.format(q,struct,q,struct))
        #break


