import os
import sys
from grouper import grouper

queue = 'rondror'
group_size = 10
cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py {} {} {} {} {}\n'

def mcss(lm, h=None):
    os.system('mkdir -p mcss')
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))
    os.system('rm -f mcss/{}/*.sh'.format(lm.sp['mcss']))

    for fname in os.listdir('mcss/{}'.format(lm.sp['mcss'])):
        if fname[-3:] == 'out':
            with open('mcss/{}/{}'.format(lm.sp['mcss'], fname)) as f:
                for line in f:
                    line = line.strip().split(' ')
                    if line[0] == 'error?':
                        print line[1]
                        #print f.read()
                        #os.system('rm mcss/{}/{}/{}*'.format(out_dir, st, line[1]))

    all_mcss = sorted([m.split('.')[0] for m in os.listdir('ligands/mcss/{}'.format(lm.sp['mcss'])) if m[-3:] == 'txt'])
    all_mcss = [tuple(m.split('-')) for m in all_mcss if m[0] != '.']
    docked = set(lm.docked(lm.pdb+lm.chembl()))

    unfinished_pairs = []
    for l1, l2 in all_mcss:
        if l1 not in docked or l2 not in docked: continue
        if os.path.exists('mcss/{}/{}-{}-to-{}.csv'.format(lm.sp['mcss'], l1, l2, lm.st)): continue
        
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
    
    os.chdir('mcss/{}'.format(lm.sp['mcss']))
    os.system('rm -f *.sh')
    for i, pairs in enumerate(grouper(group_size, unfinished_pairs)):
        with open('{}_in.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                f.write(cmd.format(lm.sp['code'], l1, l2, lm.st, lm.sp['docking'], lm.sp['mcss']))
            f.write('wait\n')
        os.system('sbatch --tasks=1 --cpus-per-task=1 -p {} -t 2:00:00 {}_in.sh'.format(queue,i))

    os.chdir('../..')




