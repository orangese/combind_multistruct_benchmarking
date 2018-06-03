import os
import sys

import itertools

out_dir = 'mcss7'
mcss_dir = 'ligands/mcss/{}'.format(out_dir)

script = '/scratch/PI/rondror/jbelk/method/combind/2_ifp/proc_mcss.py'

queue = 'rondror'

def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

def mcss(lm, h=None):
    if not os.path.exists(lm.gdir): return

    os.system('mkdir -p mcss')
    os.system('mkdir -p mcss/{}'.format(out_dir))
    os.system('rm -f mcss/{}/*.sh'.format(out_dir))

    for fname in os.listdir('mcss/{}'.format(out_dir)):
        if fname[-3:] == 'out':
            with open('mcss/{}/{}'.format(out_dir, fname)) as f:
                for line in f:
                    line = line.strip().split(' ')
                    if line[0] == 'error?':
                        print line[1]
                        #print f.read()
                        #os.system('rm mcss/{}/{}/{}*'.format(out_dir, st, line[1]))

    all_mcss = sorted([m.split('.')[0] for m in os.listdir(mcss_dir) if m.split('.')[-1] == 'txt'])
    all_mcss = [tuple(m.split('-')) for m in all_mcss]
    docked = set(lm.docked(lm.pdb+lm.chembl()))

    unfinished_pairs = []
    for l1, l2 in all_mcss:
        if l1 not in docked or l2 not in docked: continue
        if os.path.exists('mcss/{}/{}-{}-to-{}.csv'.format(out_dir, l1, l2, lm.st)): continue
        
        if h is not None:
            for f, f_data in h.items():
                for q, hlist in f_data.items():
                    hlist = set(hlist)
                    if l1 in hlist and l2 in hlist:
                        unfinished_pairs.append((l1, l2))
                        break
                if len(unfinished_pairs) == 0: continue
                if unfinished_pairs[-1] == (l1, l2): break

        if l1[:6] != 'CHEMBL' and l2[:6] != 'CHEMBL':
            unfinished_pairs.append((l1,l2))

    if len(unfinished_pairs) > 0:
        print len(unfinished_pairs), 'mcss rmsd left'
    
    group_size = 12
    os.chdir('mcss/{}'.format(out_dir))
    os.system('rm -f *.sh')
    for i, pairs in enumerate(grouper(group_size, unfinished_pairs)):
        with open('{}_in.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                f.write('$SCHRODINGER/run {} {} {} {} {} {}\n'.format(script, l1, l2, lm.st, lm.gdir, out_dir))
            f.write('wait\n')
        os.system('sbatch --tasks=6 --cpus-per-task=1 -p {} -t 1:30:00 {}_in.sh'.format(queue,i))

    os.chdir('../..')




