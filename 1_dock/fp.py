import os
import sys

import time
import itertools

FUZZY_SCRIPT = "/scratch/PI/rondror/jbelk/method/combind/2_ifp/fuzzyifp.py"

glide_dir = 'docking/glide12'
output_dir = 'ifp/ifp2'

queue = 'rondror'

def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

def get_fp(fp_list):
    if len(fp_list) > 0:
        print len(fp_list), 'fp left'
    group_size = 1
    os.chdir(output_dir)
    for i, pairs in enumerate(grouper(group_size, fp_list)):
        with open('{}fp.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                input_file = '../../{}/{}/{}_pv.maegz'.format(glide_dir, p, p)
                output_file = '{}.fp'.format(p)
                f.write('$SCHRODINGER/run {} -mode pv -input_file {} -output_file {}\n'.format(FUZZY_SCRIPT, input_file, output_file)) 
            f.write('wait\n')
        os.system('sbatch --tasks={} --cpus-per-task=1 --time=02:00:00 -p {} {}fp.sh'.format(group_size, queue, i))
    os.chdir('../..')

def structure_fp():
    os.chdir(output_dir)
    for lig in sorted(os.listdir('../../structures/ligands')):
        pdb = lig.split('_')[0]
        output_file = '{}_struct.fp'.format(pdb)
        if os.path.exists(output_file): continue
        print 'structure fp', pdb
        with open('{}.sh'.format(pdb), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('$SCHRODINGER/run {} -mode st -output_file {}\n'.format(FUZZY_SCRIPT, output_file)) 
        os.system('sbatch --time=00:10:00 -n 1 -p {} {}.sh'.format(queue, pdb))
    os.chdir('../..')

def fp(lm):
    os.system('mkdir -p ifp')
    os.system('mkdir -p {}'.format(output_dir))

    #for fp in os.listdir(output_dir):
    #    if fp.split('_')[-1] == 'struct.fp': continue
    #    if not os.path.exists('{}/{}'.format(glide_dir, fp.split('.')[0])):
            #os.system('rm -f {}/{}'.format(output_dir,fp))
            #print 'glide folder has been deleted since fp creation'
            #print 'removing fp', fp

    unfinished = []    
    for lig in lm.docked(lm.pdb+lm.chembl()):
        output_file = '{}/{}-to-{}.fp'.format(output_dir, lig, lm.st)
        if os.path.exists(output_file): continue
        unfinished.append('{}-to-{}'.format(lig, lm.st))
        
    structure_fp()
    get_fp(unfinished)



