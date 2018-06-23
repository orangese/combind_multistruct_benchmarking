import os
import sys
from grouper import grouper

from schrodinger.structure import StructureReader, StructureWriter

queue = 'rondror'

rmsd_cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py RMSD {} {} {}\n'
size_cmd = '$SCHRODINGER/run {}/2_ifp/mcss_main.py SIZE {}\n'
init_cmd = '$SCHRODINGER/utilities/canvasMCS -imae {}-{}_in.mae -ocsv {} -stop 10 {} -atomtype C {}/2_ifp/custom_types/{}.typ\n'

opt = ['']#,'-nobreakaring']

class SM:
    def __init__(self, lm):
        self.mdir, self.gdir, self.st = lm.sp['mcss'], lm.sp['docking'], lm.st
        self.tf = lm.sp['mcss_type']
        self.lm = lm
        self.docked = set(lm.docked(lm.pdb+lm.chembl()))
        self.no_mcss = set([])
        self.no_size = set([])
        self.no_rmsd = set([])
    
    def get_path(self,l1,l2,o,add_dir=False,add_st=False,ext='csv'):
        pth = '{}-{}-{}{}'.format(l1,l2,self.tf,o)
        if add_st: pth = '{}-{}-{}'.format(pth,self.st,self.gdir)
        if add_dir: pth = 'mcss/{}/{}-{}/{}'.format(self.mdir,l1,l2,pth)
        return '{}.{}'.format(pth,ext)

    def add(self, l1, l2, rmsd=False):
        for o in opt:
            if not os.path.exists(self.get_path(l1,l2,o,add_dir=True)):
                self.no_mcss.add((l1,l2,o))
            elif not os.path.exists(self.get_path(l1,l2,o,add_dir=True,ext='size')):
                self.no_size.add((l1,l2,o))
            elif rmsd and l1 in self.docked and l2 in self.docked:
                if not os.path.exists(self.get_path(l1,l2,o,add_dir=True,add_st=True)):
                    self.no_rmsd.add((l1,l2,o))

        if not os.path.exists('mcss/{}/{}-{}/{}-{}_in.mae'.format(self.mdir,l1,l2,l1,l2)):
            os.system('mkdir -p mcss/{}/{}-{}'.format(self.mdir,l1,l2))
            stwr = StructureWriter('mcss/{}/{}-{}/{}-{}_in.mae'.format(self.mdir,l1,l2,l1,l2))
            stwr.append(StructureReader('ligands/prepared_ligands/{}/{}.mae'.format(l1,l1)).next())
            stwr.append(StructureReader('ligands/prepared_ligands/{}/{}.mae'.format(l2,l2)).next())
            stwr.close()

    def proc(self, init=False, rmsd=False, size=False):
        if init: 
            group_size = 100
            all_pairs = self.no_mcss
        if size:
            group_size = 25
            all_pairs = self.no_size
        if rmsd:
            group_size = 10
            all_pairs = self.no_rmsd

        os.chdir('mcss/{}'.format(self.mdir))
        for i, pairs in enumerate(grouper(group_size, all_pairs)):
            if init: script = 'init{}.sh'.format(i)
            if size: script = 'size{}.sh'.format(i)
            if rmsd: script = 'rmsd{}.sh'.format(i)
            with open(script, 'w') as f:
                f.write('#!/bin/bash\nmodule load schrodinger\n')
                if init: f.write('export SCHRODINGER_CANVAS_MAX_MEM=1e+12\n')
                for p in pairs:
                    if p is None: continue
                    l1,l2,o = p
                    f.write('cd {}-{}\n'.format(l1,l2))
                    if init: f.write(init_cmd.format(l1,l2,self.get_path(l1,l2,o),o,self.lm.sp['code'], self.tf))
                    if size: f.write(size_cmd.format(self.lm.sp['code'], self.get_path(l1,l2,o)))
                    #f.write(init_cmd.format(l1,l2,self.get_path(l1,l2,o),o,self.lm.sp['code'], self.tf))
                    #f.write(size_cmd.format(self.lm.sp['code'], self.get_path(l1,l2,o)))
                    if rmsd: f.write(rmsd_cmd.format(self.lm.sp['code'], self.get_path(l1,l2,o),self.st,self.gdir))
                    f.write('cd ..\n')
                f.write('wait\n')
            os.system('sbatch -p {} --tasks=1 --cpus-per-task=1 --nice -t 8:00:00 {}'.format(queue,script))
        os.chdir('../..')


def mcss(lm, chembl={}, max_num=20):
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))
    
    sm = SM(lm)

    for f, f_data in chembl.items():
        for q,c in f_data.items():
            #print q,c
            for i,l1 in enumerate(c):
                if l1 == '': continue
                sm.add(q,l1,True)
                for l2 in c[i+1:]: #pass
                    if l1 < l2: sm.add(l1,l2,True)
                    else: sm.add(l2,l1,True)

    for i,l1 in enumerate(lm.pdb[:max_num]):
        for l2 in lm.pdb[i+1:max_num]:
            sm.add(l1,l2,True)

    for l1 in lm.pdb[:max_num]:
        for l2 in lm.chembl():
            sm.add(l1,l2)
    
    if len(sm.no_mcss) > 0:
        print len(sm.no_mcss), 'mcss init pairs left'
        sm.proc(init=True)
    if len(sm.no_size) > 0:
        print len(sm.no_size), 'mcss size pairs left'
        sm.proc(size=True)
    if len(sm.no_rmsd) > 0:
        print len(sm.no_rmsd), 'mcss rmsd pairs left'
        sm.proc(rmsd=True)





