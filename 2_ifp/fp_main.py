import os
import sys
import time

FUZZY_SCRIPT = "/scratch/PI/rondror/jbelk/method/code/2_ifp/fuzzyifp.py"

glide_dir = sys.argv[1]
output_dir = sys.argv[2]

d_list = sys.argv[3:]

os.chdir('../../data')

if len(d_list) == 0:
    d_list = sorted(os.listdir('.'))

queue = 'rondror'

#datasets = {'B1AR':'2VT4','A2AR':'2YDO'}

datasets = {'D2R':'6CM4','AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
            'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23'}

for d in sorted(d_list):
    #if d[0] == '.': continue
    #if d == 'DOR_allwat': continue

    if d not in datasets: continue
    st = datasets[d]

    print d
    
    if not os.path.exists('{}/docking/{}'.format(d, glide_dir)): continue

    #os.system('rm -rf {}/ifp'.format(d))
    os.system('mkdir -p {}/ifp'.format(d))
    os.chdir('{}/ifp'.format(d))
    os.system('mkdir -p ' + output_dir)
    os.chdir(output_dir)

    os.system('rm -f *.out *.sh')
    for fp in os.listdir('.'):
        if os.path.getsize(fp) == 0:
            print 'found empty file', fp
            os.system('rm -f {}'.format(fp))
        if len(fp) == 7: continue
        if not os.path.exists('../../docking/{}/{}'.format(glide_dir, fp.split('.')[0])):
            os.system('rm -f {}'.format(fp))
            print 'glide folder has been deleted since fp creation'
            print 'removing fp', fp

    #os.chdir('../../..')
    #continue
    for lig in sorted(os.listdir('../../structures/ligands')):
        pdb = lig.split('_')[0]
        output_file = '{}.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print pdb
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('$SCHRODINGER/run {} -mode st -output_file {}\n'.format(FUZZY_SCRIPT, output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))
    
    dock_path = '../../docking/{}'.format(glide_dir)
    for dock in sorted(os.listdir(dock_path)):
        input_file = '../../docking/{}/{}/{}_pv.maegz'.format(glide_dir, dock, dock)
        output_file = '{}.fp'.format(dock)
        if not os.path.exists(input_file): continue
        if os.path.exists(output_file): continue
        l, s = dock.split('-to-')
        if s != st: continue
        print dock
        with open('{}.sh'.format(dock), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('$SCHRODINGER/run {} -mode pv -input_file {} -output_file {}\n'.format(FUZZY_SCRIPT, input_file, output_file)) 
        os.system('sbatch --time=03:00:00 -n 1 -p {} {}.sh'.format(queue, dock))

    os.chdir('../../..')

