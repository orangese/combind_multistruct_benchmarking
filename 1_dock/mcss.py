import os
import sys
from grouper import grouper

csv_out  = 'mcss/{}/{}/{}.csv'
mae_out  = 'mcss/{}/{}/{}_mcss.mae'
rmsd_out = 'mcss/{}/{}/{}-to-{}-{}.csv'

queue = 'rondror'

type_file = 'mcss1'

init_cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py INIT {} {} {}\n'
rmsd_cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py RMSD {} {} {} {}\n'

def mcss(lm, chembl=None, max_num=20):
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))

    all_pairs = set([])
    if chembl is not None:
        for f, f_data in chembl.items():
            for q,c in f_data.items():
                for i in range(len(c)):
                    for j in range(i+1,len(c)):
                        all_pairs.add((c[i],c[j]))

    pdb_ligs = lm.pdb[:max_num]
    for i1 in range(len(pdb_ligs)):
        for i2 in range(i1+1,len(pdb_ligs)):
            l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
            all_pairs.add((l1, l2))

    for l1 in pdb_ligs:
        for l2 in lm.chembl():
            all_pairs.add((l1, l2))

    mdir = lm.sp['mcss']

    no_mcss = []
    no_rmsd = []
    docked = set(lm.docked(lm.pdb+lm.chembl()))
    for l1,l2 in all_pairs:
        name = '{}-{}'.format(l1,l2)
        if os.path.exists(rmsd_out.format(mdir,name,name,lm.st,lm.sp['docking'])): 
            continue
        if not os.path.exists(csv_out.format(mdir,name,name)) or not os.path.exists(mae_out.format(mdir,name,name)):
            os.system('rm -rf mcss/{}/{}'.format(mdir,name))
            no_mcss.append((l1,l2))
        elif l1 in docked and l2 in docked:
            no_rmsd.append((l1,l2))

    if len(no_mcss) > 0:
        print len(no_mcss), 'mcss init pairs left'
        proc(no_mcss, lm, init=True, rmsd=False)
    if len(no_rmsd) > 0:
        print len(no_rmsd), 'mcss rmsd pairs left'
        proc(no_rmsd, lm, init=False, rmsd=True)

def proc(all_pairs, lm, init=False, rmsd=False):
    if init: group_size = 100
    if rmsd: group_size = 10
    for i, pairs in enumerate(grouper(group_size, all_pairs)):
        if init: script = 'init{}.sh'.format(i)
        if rmsd: script = 'rmsd{}.sh'.format(i)
        with open('mcss/{}/{}'.format(lm.sp['mcss'],script), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2)
                if init: f.write('mkdir -p {}\n'.format(name))
                f.write('cd {}\n'.format(name))
                if init: f.write(init_cmd.format(lm.sp['code'], l1, l2, type_file))
                if rmsd: f.write(rmsd_cmd.format(lm.sp['code'], l1, l2, lm.st, lm.sp['docking']))
                f.write('cd ..\n')
            f.write('wait\n')
        os.chdir('mcss/{}'.format(lm.sp['mcss']))
        os.system('sbatch -p {} --tasks=1 --cpus-per-task=1 -t 2:00:00 {}'.format(queue,script))
        os.chdir('../..')





