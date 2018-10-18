import os
import sys
from shared_paths import shared_paths
from grouper import grouper

queue = 'owners'
group_size=10

# TODO: features seperately set in statistics.py
features = ['sb1','sb2','sb3','mcss','hbond','pipi','contact']

def compute(lm, all_pairs):
    os.system('mkdir -p stats')
    os.system('mkdir -p stats/{}'.format(shared_paths['stats']))
    os.chdir('stats/{}'.format(shared_paths['stats']))
    
    for i,group in enumerate(grouper(group_size, all_pairs)):
        with open('stats{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            f.write('unset PYTHONPATH\n')
            f.write('unset PYTHONHOME\n')
            for pair in group:
                if pair is None: continue
                l1,l2 = pair
                f.write('python {}/3_analyze/statistics.py {} {} {}\n'.format(shared_paths['code'],
                    lm.prot, l1, l2))
        os.system('sbatch -p {} -t 1:00:00 stats{}.sh'.format(queue, i))
    os.chdir('../..')

def stats(lm, max_ligs = 20):
    """
    Compute statistics for lm.prot.
    """
    ligs = lm.docked(lm.pdb)
    num_ligs = min(max_ligs, len(ligs))
    not_done = set()

    for i in range(num_ligs):
        for j in range(i+1,num_ligs):
            l1,l2 = ligs[i],ligs[j]
            for k in features:
                if not os.path.exists('stats/{}/{}-{}-to-{}-{}.txt'.format(shared_paths['stats'],l1,l2,lm.st,k)):
                    print(l1,l2,k)
                    not_done.add((l1,l2))

    if len(not_done) > 0:
        print(len(not_done), 'stats pairs left')
        compute(lm, not_done)
