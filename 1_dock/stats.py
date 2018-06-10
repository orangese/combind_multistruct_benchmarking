import os
import sys

import itertools

group_size=5
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

queue = 'rondror'
out_dir = 'stats/stats2'
features = ['sb1','sb2','sb3','mcss','hbond','pipi','contact']

stats_script = '/scratch/PI/rondror/jbelk/method/combind/3_analyze/statistics.py'

def stats(lm):
    os.system('rm -rf {}'.format(out_dir))
    ligs = lm.docked(lm.pdb)
    num_ligs = min(10, len(ligs))
    #print ligs[:num_ligs]
    not_done = []
    for i in range(num_ligs):
        for j in range(i+1,num_ligs):
            l1,l2 = ligs[i],ligs[j]
            for k in features:
                if not os.path.exists('{}/{}-{}-to-{}-{}.txt'.format(out_dir,l1,l2,lm.st,k)):
                    not_done.append((l1,l2))
                    break

    if len(not_done) > 0:
        print len(not_done), 'stats pairs left'
        compute(lm, not_done)

def compute(lm, all_pairs):
    #print all_pairs
    os.system('mkdir -p stats')
    os.system('mkdir -p {}'.format(out_dir))
    os.chdir(out_dir)
    
    for i,group in enumerate(grouper(group_size, all_pairs)):
        with open('stats{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            for pair in group:
                if pair is None: continue
                l1,l2 = pair
                f.write('python {} {} {} {}\n'.format(stats_script, lm.prot, l1, l2))
            f.write('rm stats{}.sh\n'.format(i))
        os.system('sbatch -p {} -t 1:00:00 stats{}.sh'.format(queue, i))

    os.chdir('../..')

