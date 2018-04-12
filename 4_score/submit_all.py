import os
import sys

sys.path.append('../3_analyze')
from containers import Dataset

output_dir = 'scores/test_pdb2'
script = '/scratch/PI/rondror/jbelk/method/code/4_score/score_query.py'

struct_dict = {'AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
               'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23'}

settings = {
    'data_dir' : '/scratch/PI/rondror/jbelk/method/data',
    'glide_dir' : 'docking/glide12',
    'ifp_dir' : 'ifp/ifp3',
    'mcss_dir' : 'mcss/mcss1',
    'features' : {'mcss':[],'hbond':[2,3]},#,'pipi':[6],'contact':[11]},
    'stats_prots' : ['CHK1','B1AR','ERA'],
    'struct_dict' : struct_dict,
    'num_stats_ligs' : 10,
    'num_stats_poses' : 25,
    'smooth' : 0.5,
    'num_pred_chembl' : 12,
    'num_poses' : 100,
    't' : 10,
    'req_stereo':True,
    'use_chembl':False
}

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def write_job_file(out_path, query, prot):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('python {} {} {}\n'.format(script, query, prot))

predict_prots = ['BRAF']
predict_data = Dataset(predict_prots, struct_dict, settings['data_dir'], 
                       settings['glide_dir'], settings['ifp_dir'], settings['mcss_dir'])

for p, prot in predict_data.proteins.items():
    os.system('mkdir -p {}/{}/{}'.format(settings['data_dir'], p, output_dir))
    os.chdir('{}/{}/{}'.format(settings['data_dir'], p, output_dir))
    write_settings_file('settings.py', settings)
    for q in prot.lm.get_pdb():
        write_job_file('{}.in'.format(q), q, p)
        os.system('sbatch -p owners -t 1:00:00 -o {}.out {}.in'.format(q,q))
        #break
    #break


