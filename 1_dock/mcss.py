import os
import sys
from grouper import grouper

from schrodinger.structure import StructureReader, StructureWriter

csv_out  = 'mcss/{}/{}/{}.csv'
rmsd_out = 'mcss/{}/{}/{}-to-{}-{}.csv'

queue = 'rondror'

rmsd_cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py RMSD {} {} {} {}\n'
init_cmd = '$SCHRODINGER/utilities/canvasMCS -imae {}_in.mae -ocsv {}.csv -stop 10 -atomtype C {}/2_ifp/custom_types/mcss1.typ\n'

class SM:
    def __init__(self, lm):
        self.mdir, self.gdir, self.st = lm.sp['mcss'], lm.sp['docking'], lm.st
        self.docked = set(lm.docked(lm.pdb+lm.chembl()))
        self.no_mcss = set([])
        self.no_rmsd = set([])

    def rmsd_done(self, pair):
        return os.path.exists(rmsd_out.format(self.mdir, pair, pair, self.st, self.gdir))

    def init_done(self, pair):
        return os.path.exists(csv_out.format(self.mdir,pair,pair))

    def add(self, l1, l2, rmsd=False):
        if not self.init_done('{}-{}'.format(l1,l2)):
            self.no_mcss.add((l1,l2))
        elif rmsd and l1 in self.docked and l2 in self.docked and not self.rmsd_done('{}-{}'.format(l1,l2)):
            self.no_rmsd.add((l1,l2))

def mcss(lm, chembl={}, max_num=20):
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))
    
    sm = SM(lm)

    for f, f_data in chembl.items():
        for q,c in f_data.items():
            for i,l1 in enumerate(c):
                for l2 in c[i+1:]:
                    sm.add(l1,l2,True)

    for i,l1 in enumerate(lm.pdb[:max_num]):
        for l2 in lm.pdb[i+1:max_num]:
            sm.add(l1,l2,True)

    for l1 in lm.pdb[:max_num]:
        for l2 in lm.chembl():
            sm.add(l1,l2)

    if len(sm.no_mcss) > 0:
        print len(sm.no_mcss), 'mcss init pairs left'
        proc(sm.no_mcss, lm, init=True, rmsd=False)
    if len(sm.no_rmsd) > 0:
        print len(sm.no_rmsd), 'mcss rmsd pairs left'
        proc(sm.no_rmsd, lm, init=False, rmsd=True)

def proc(all_pairs, lm, init=False, rmsd=False):
    if init: group_size = 100
    if rmsd: group_size = 5
    for i, pairs in enumerate(grouper(group_size, all_pairs)):
        if init: script = 'init{}.sh'.format(i)
        if rmsd: script = 'rmsd{}.sh'.format(i)
        with open('mcss/{}/{}'.format(lm.sp['mcss'],script), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2)
                f.write('cd {}\n'.format(name))
                if init:
                    os.system('rm -rf mcss/{}/{}'.format(lm.sp['mcss'],name))
                    os.system('mkdir mcss/{}/{}'.format(lm.sp['mcss'],name))
                    stwr = StructureWriter('mcss/{}/{}/{}_in.mae'.format(lm.sp['mcss'],name,name))
                    stwr.append(StructureReader('ligands/prepared_ligands/{}/{}.mae'.format(l1,l1)).next())
                    stwr.append(StructureReader('ligands/prepared_ligands/{}/{}.mae'.format(l2,l2)).next())
                    stwr.close()
                    f.write(init_cmd.format(name, name, lm.sp['code']))
                    f.write('rm -f {}_in.mae\n'.format(name))
                if rmsd: f.write(rmsd_cmd.format(lm.sp['code'], l1, l2, lm.st, lm.sp['docking']))
                f.write('cd ..\n')
            f.write('wait\n')
        os.chdir('mcss/{}'.format(lm.sp['mcss']))
        os.system('sbatch -p {} --tasks=1 --cpus-per-task=1 -t 2:00:00 {}'.format(queue,script))
        os.chdir('../..')





