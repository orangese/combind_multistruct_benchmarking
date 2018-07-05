import os
import sys
from grouper import grouper

queue = 'owners'
group_size = 10

pv_cmd = '$SCHRODINGER/run {}/2_ifp/fp_main.py -mode pv -input_file {} -output_file {}\n'
st_cmd = '$SCHRODINGER/run {}/2_ifp/fp_main.py -mode st -output_file {}\n'

def get_fp(lm, fp_list):
    if len(fp_list) > 0:
        print len(fp_list), 'fp left'
    for i, pairs in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                input_file = '../../docking/{}/{}/{}_pv.maegz'.format(lm.sp['docking'], p, p)
                output_file = '{}.fp'.format(p)
                f.write(pv_cmd.format(lm.sp['code'], input_file, output_file)) 
            f.write('wait\n')
        os.system('sbatch --cpus-per-task=1 --time=02:00:00 -p {} {}fp.sh'.format(queue, i))

def structure_fp(lm):
    for lig in sorted(os.listdir('../../structures/ligands')):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print 'structure fp', pdb
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write(st_cmd.format(lm.sp['code'], output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))

def fp(lm):
    os.system('mkdir -p ifp')
    os.system('mkdir -p ifp/{}'.format(lm.sp['ifp']))

    unfinished = []    
    for lig in lm.docked(lm.pdb+lm.chembl()):
        output_file = 'ifp/{}/{}-to-{}.fp'.format(lm.sp['ifp'], lig, lm.st)
        if os.path.exists(output_file): continue
        unfinished.append('{}-to-{}'.format(lig, lm.st))
       
    os.chdir('ifp/{}'.format(lm.sp['ifp'])) 
    structure_fp(lm)
    get_fp(lm, unfinished)
    os.chdir('../..')
