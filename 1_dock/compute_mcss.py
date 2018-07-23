"""
Controls computation of MCSS features.
"""

import os
import sys
from grouper import grouper
from schrodinger.structure import StructureReader, StructureWriter

class MCSSController:
    """
    Records Ligand pairs with incomplete MCSS features and schedules
    their computation.

    lm: LigandManager
    max_pdb: int, maximum number of pdb ligands to consider
    docked: set(string), Set of ligands with docking results
    no_mcss: set(string), Ligands with no init MCSS
    no_size: set(string), Ligands with no size MCSS
    no_rmsd: set(string), Ligands with no MCSS RMSD
    """

    INIT_GROUP_SIZE = 100
    SIZE_GROUP_SIZE = 25
    RMSD_GROUP_SIZE = 10

    QUEUE = 'owners'
    
    TEMPLATE = ('#!/bin/bash\n'
                '#SBATCH -p {}\n'
                '#SBATCH --tasks=1\n'
                '#SBATCH -t 8:00:00\n'
                'ml load chemistry\n'
                'ml load schrodinger/2018-1\n'
                '{}'
                'wait\n'
                ).format(QUEUE, '{}')

    def __init__(self, lm, max_pdb, max_poses):
        self.lm = lm
        self.max_pdb = max_pdb
        
        self.pdb     = self.lm.docked(self.lm.pdb)[:self.max_pdb]
        self.docked  = set(lm.docked(lm.pdb+lm.chembl()))
        self.no_mcss = set([])
        self.no_size = set([])
        self.no_rmsd = set([])

        # Remaining '{}' to be filled with mode, ligand1, ligand2
        self.command = ('$SCHRODINGER/run {0:}/2_ifp/mcss_main.py {1:} '
                        '{1:} {1:} {2:} {3:} {4:} {5:} {6:} {7:} {8:}'
                        ).format(lm.sp['code'], '{}', max_poses, lm.st, lm.prot,
                        lm.sp['mcss'], lm.sp['docking'], lm.sp['data'], lm.sp['code'])

    # Methods to add ligand pairs
    def add_pdb_to_pdb(self):
        """
        Add all pdb to pdb ligand pairs
        """
        for i, l1 in enumerate(self.pdb):
            for l2 in self.pdb[i+1:]:
                self._add(l1, l2, True)

    def add_pdb_to_allchembl(self):
        """
        Add all pdb - chembl ligand pairs.
        Most of these are only used for deciding which chembl ligands to use
        so don't compute rmsd for all of them.
        """
        for l1 in self.pdb:
            for l2 in self.lm.chembl():
                self._add(l1,l2, False)

    def add_pick_helpers(self, pick_helpers):
        """
        Adds pairs of pdb to chembl ligands and chembl to chembl ligands
        that will be jointly used in a score computation as specified
        by pick_helpers.

        Chembl - Chembl pairs that are specified by "chembl"
        pick_helpers: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
        """
        for fname, queries in pick_helpers.items():
            for query, chembl_ligands in queries.items():
                for i, l1 in enumerate(chembl_ligands):
                    assert l1 != '' # switch to continue if this happens
                    self._add(query, l1, True)
                    for l2 in chembl_ligands[i+1:]:
                        assert l2 != '' # switch to continue if this happens
                        self._add(l1, l2, True)

    def _add(self, l1, l2, rmsd):
        """
        l1: string, ligand name
        l2: string, ligand name
        rmsd: bool, compute RMSD if true
        """
        if l1 > l2: l1, l2 = l2, l1

        if not os.path.exists(self._path_to_mcss(l1, l2)):
            os.system('mkdir -p {}'.format(self._path_to_directory(l1, l2)))
            self.no_mcss.add((l1,l2))
        elif not os.path.exists(self._path_to_mcss(l1, l2, ext='size')):
            self.no_size.add((l1,l2))
        elif (rmsd and l1 in self.docked and l2 in self.docked
              and not os.path.exists(self._path_to_mcss(l1, l2, add_docking=True))):
            self.no_rmsd.add((l1, l2))

    # Methods to execute computation
    def execute(self):
        """
        Execute all jobs.
        """
        if len(self.no_mcss) > 0:
            print(len(self.no_mcss), 'mcss init pairs left')
            self._execute('INIT')
        if len(self.no_size) > 0:
            print(len(self.no_size), 'mcss size pairs left')
            self._execute('SIZE')
        if len(self.no_rmsd) > 0:
            print(len(self.no_rmsd), 'mcss rmsd pairs left')
            self._execute('RMSD')

    def _execute(self, mode):
        """
        Schedules execution of all jobs with mode "mode".

        mode: string in ['INIT', 'SIZE', 'RMSD']
        """

        if mode == 'INIT': 
            group_size = self.INIT_GROUP_SIZE
            all_pairs = self.no_mcss
            header = 'export SCHRODINGER_CANVAS_MAX_MEM=1e+12\n'
            script = 'init{}.sh'
        elif mode == 'SIZE':
            group_size = self.SIZE_GROUP_SIZE
            all_pairs = self.no_size
            header = ''
            script = 'size{}.sh'
        elif mode == 'RMSD':
            group_size = self.RMSD_GROUP_SIZE
            all_pairs = self.no_rmsd
            header = ''
            script = 'rmsd{}.sh'
        else:
            assert False, "Mode {} not in ['INIT', 'SIZE', 'RMSD']".format(mode)

        os.chdir('mcss/{}'.format(self.lm.sp['mcss']))
        for i, pairs in enumerate(grouper(group_size, all_pairs)):
            contents = header
            for l1, l2 in pairs:
                contents += '( cd {}-{}; '.format(l1,l2)
                contents += self.command.format(mode, l1, l2)
                contents += ' ; cd .. ) \n'
            with open(script.format(i), 'w') as f: f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script.format(i)))
        os.chdir('../..')

    # Methods to get file paths
    def _name(self, l1, l2):
        """
        l1: string, ligand name
        l2: string, ligand name
        """
        if l1 > l2: l1, l2 = l2, l1
        return "{}-{}".format(l1, l2)

    def _path_to_directory(self, l1, l2):
        """
        l1: string, ligand name
        l2: string, ligand name
        """
        return 'mcss/{}/{}'.format(self.lm.sp['mcss'], self._name(l1, l2))

    def _path_to_mcss(self, l1, l2, add_docking=False, ext='csv'):
        """
        l1: string, ligand name
        l2: string, ligand name

        add_dir: bool, whether to add path from dataset root
        add_docking: bool, whether to specify docking run
        """

        path = self._name(l1, l2)

        if add_docking: path = '{}-{}-{}'.format(path, self.lm.st, self.lm.sp['docking'])
        path = '{}/{}'.format(self._path_to_directory(l1, l2), path)
        path = '{}.{}'.format(path, ext)
        return path

def mcss(lm, pick_helpers={}, max_pdb=20, max_poses=100):
    """
    Compute unfinished MCSS features.

    lm: LigandManager instance
    chembl: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
    max_pdb: int, maximum number of pdb ligands to consider
    """
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))
    controller = MCSSController(lm, max_pdb, max_poses)
    
    controller.add_pdb_to_pdb()
    controller.add_pdb_to_allchembl()
    controller.add_pick_helpers(pick_helpers)
    
    controller.execute()
