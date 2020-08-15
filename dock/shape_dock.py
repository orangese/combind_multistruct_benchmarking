import os
import sys
from utils import grouper

QUEUE = 'owners'
GROUP_SIZE = 25

QUERY = '../../ligands/{query}/{query}.mae'
REF = '../../structures/ligands/{ref}_lig.mae'
NATIVE = '../../structures/ligands/{query}.mae'
POSE = '{query}-to-{ref}/{ref}_lig_align.maegz'
RMSD = '{query}-to-{ref}/rmsd.csv'

dock_cmd = '$SCHRODINGER/shape_screen -shape {REF} -screen {QUERY} -WAIT -flex -atomtypes mmod -report 1000 -sample thorough\n'
rmsd_cmd = 'run rmsd.py -use_neutral_scaffold -pv second -c rmsd.csv {NATIVE} {POSE}\n'

def get_state(ref, query):
    # 0: do nothing
    # 1: compute rmsd
    # 2: dock
    if os.path.exists(RMSD.format(ref=ref, query=query)):
        return 0
    
    if os.path.exists(POSE.format(ref=ref, query=query)):
        if os.path.exists(NATIVE.format(query=query)): return 1
        else: return 0
    
    if not (    os.path.exists(QUERY.format(query=query))
            and os.path.exists(REF.format(ref=ref))):
        print(REF.format(ref=ref))
        print(QUERY.format(query=query))
        return 0
    return 2

def proc_all(all_pairs, dock=False, rmsd=False):
    os.system('rm sdock*.sh')
    for i,group in enumerate(grouper(GROUP_SIZE, all_pairs)):
        with open('sdock{}.sh'.format(i),'w') as f:
            f.write('#!/bin/sh\n')
            for pair in group:
                if pair is None: continue
                ref, query = pair

                f.write('cd {}-to-{}\n'.format(query, ref))
                if dock:
                    os.system('rm -rf {}-to-{}'.format(query, ref))
                    os.system('mkdir {}-to-{}'.format(query, ref))
                    f.write(dock_cmd.format(REF='../'+REF.format(ref=ref),
                                            QUERY='../'+QUERY.format(query=query)
                                            ))
                if rmsd:
                    f.write(rmsd_cmd.format(NATIVE='../'+NATIVE.format(query=query),
                                            POSE='../'+POSE.format(ref=ref, query=query)))
                f.write('cd ..\n')

            f.write('which shape_screen')

        os.system('sbatch -p {} -t 1:00:00 -o sdock.out sdock{}.sh'.format(QUEUE, i))

def shape_dock(lm):
    if lm.st is None:
        return

    docking = 'shape'
    ref = lm.st
    queries = lm.get_pdb()

    os.system('mkdir -p docking/{}'.format(docking))
    os.chdir('docking/{}'.format(docking))

    to_dock = []
    to_rmsd = []
    for query in queries:
        s = get_state(ref, query)
        if s == 2: to_dock.append((ref, query))
        if s == 1: to_rmsd.append((ref, query))

    if len(to_dock) > 0:
        print('shape docking {} ligands'.format(len(to_dock)))
        proc_all(to_dock, dock=True)
    if len(to_rmsd) > 0:
        print('computing {} rmsds'.format(len(to_rmsd)))
        proc_all(to_rmsd, rmsd=True)

    os.chdir('../..')
