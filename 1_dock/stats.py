import os
import sys
from grouper import grouper

queue = 'owners'
group_size=5

features = ['sb1','sb2','sb3','mcss','hbond','pipi','contact']

def stats(lm):
    """
    Compute statistics for 10 ligands.
    """
    ligs = lm.docked(lm.pdb)
    num_ligs = min(10, len(ligs))
    not_done = []

    for i in range(num_ligs):
        for j in range(i+1,num_ligs):
            l1,l2 = ligs[i],ligs[j]
            for k in features:
                if not os.path.exists('stats/{}/{}-{}-to-{}-{}.txt'.format(lm.sp['stats'],l1,l2,lm.st,k)):
                    print l1,l2,k
                    not_done.append((l1,l2))

    if len(not_done) > 0:
        print len(not_done), 'stats pairs left'
        compute(lm, not_done)

def compute(lm, all_pairs):
    os.system('mkdir -p stats')
    os.system('mkdir -p stats/{}'.format(lm.sp['stats']))
    os.chdir('stats/{}'.format(lm.sp['stats']))
    
    for i,group in enumerate(grouper(group_size, all_pairs)):
        with open('stats{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            for pair in group:
                if pair is None: continue
                l1,l2 = pair
                f.write('python {}/3_analyze/statistics.py {} {} {}\n'.format(lm.sp['code'], lm.prot, l1, l2))
            #f.write('rm stats{}.sh\n'.format(i))
        os.system('sbatch -p {} -t 1:00:00 stats{}.sh'.format(queue, i))
    os.chdir('../..')
