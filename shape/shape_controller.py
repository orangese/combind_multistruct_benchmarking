import os
import numpy as np
from glob import glob

from utils import grouper
from shape.shape import shape

class ShapeController:
    GROUP_SIZE = 5
    QUEUE = 'owners'
    TEMPLATE = ('#!/bin/bash\n'
                '#SBATCH -p {}\n'
                '#SBATCH --tasks=1\n'
                '#SBATCH -t 2:00:00\n'
                '{}'
                ).format(QUEUE, '{}')

    def __init__(self, lm):
        self.lm = lm
        self.shapes = {}
        
        self.todo = set([])
        self.shape_file = "{}/{}-{}-{}.npy".format(self.lm.path('SHAPE'), '{0:}-{1:}', self.lm.st, lm.params['docking_version'])
        self.cmd = 'python $COMBINDHOME/shape/shape.py {2:} {3:}'
        self.cmd += ' '+self.shape_file
        self.cmd += ' --version '+self.lm.params['shape_version']
        self.cmd += ' --max-poses '+str(self.lm.params['max_poses'])

    def get(self, l1, l2, pose1, pose2):
        if l1 > l2:
            l1, l2 = l2, l1
            pose1, pose2 = pose2, pose1
        return self.shapes[(l1, l2)][pose1, pose2]

    def load(self, ligands):
        ligands = list(set(ligands))
        for i, _l1 in enumerate(ligands):
            for _l2 in ligands[i+1:]:
                if _l1 > _l2:
                    l1, l2 = _l2, _l1
                else:
                    l1, l2 = _l1, _l2
                self.shapes[(l1, l2)] = np.load(self.shape_file.format(l1, l2))

    def compute(self, ligands):
        previous_cwd = os.getcwd()
        os.system('mkdir -p {}'.format(self.lm.path('SHAPE')))
        os.chdir(self.lm.path('SHAPE'))
        if type(ligands[0]) == list:
            for _ligands in ligands:
                self._add_ligands(_ligands)
        else:
            self._add_ligands(ligands)
        self._execute()
        os.chdir(previous_cwd)

    def _add_ligands(self, ligands):
        compute_rmsds = True
        for i, l1 in enumerate(ligands):
            for l2 in ligands[i+1:]:
                self._add_pair(l1, l2)

    def _add_pair(self, l1, l2):
        if l1 > l2:
            l1, l2 = l2, l1
        fname = self.shape_file.format(l1, l2)
        docked = self.lm.docked([l1]) and self.lm.docked([l2])
        if docked and not os.path.exists(fname):
            self.todo.add((l1, l2))

    def _execute(self):
        for i, pairs in enumerate(grouper(self.GROUP_SIZE, self.todo)):
            script = 'shape{}.sh'.format(i)
            contents = ''
            for l1, l2 in pairs:
                pv1 = self.lm.path('DOCK_PV', {'ligand': l1})
                pv2 = self.lm.path('DOCK_PV', {'ligand': l2})
                contents += self.cmd.format(l1, l2, pv1, pv2)
                contents += '\n'
            with open(script, 'w') as f:
                f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))
