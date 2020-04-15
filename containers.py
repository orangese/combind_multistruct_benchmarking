"""
These classes keep track of
Poses, Ligands, Docking results for a given grid,
ligand names and preparation progress (LigandManager),
and ligands + docking results for a Protein.
"""

import os
import sys
import numpy as np
from glob import glob
import pandas as pd
from schrodinger.structure import SmilesStructure
import re
from ifp.fp_controller import parse_fp_file
from mcss.mcss_controller import MCSSController

class Pose:
    def __init__(self, rmsd, gscore, fp):
        """
        rmsd (float): RMSD to crystallographic pose
        gscore (float): glide score
        fp ({(int, string): float}): dict mapping (interactiontype, resname)
            to interaction score
        """
        self.rmsd = rmsd
        self.gscore = gscore
        self.fp = fp

class Ligand:
    def __init__(self, ligand, params, paths):
        self.ligand = ligand
        self.params = {k:v for k, v in params.items()}
        self.params['ligand'] = ligand
        self.params['pdb'] = ligand.split('_')[0]
        self.paths = paths

        self.poses = None

    def load_poses(self, load_fp):
        gscores, emodels, rmsds = self.parse_glide_output()

        fps = {}
        if load_fp:  
            fps = parse_fp_file(self.path('IFP'))

        self.poses = [Pose(rmsds[i], gscores[i], fps.get(i, {}))
                      for i in range(len(gscores))]

    def parse_glide_output(self):
        if not os.path.exists(self.path('DOCK')):
            return [], [], []

        gscores, rmsds = [], []
        with open(self.path('DOCK_REPT')) as fp:
            for line in fp:
                line = line.strip().split()
                if len(line) <= 1: continue
                lig, _lig = self.params['ligand'], line[1]
                if not (_lig == lig or _lig == lig+'_out' or _lig == '1'):
                    continue

                if _lig == '1':
                    score = line[2]
                else:
                    score = line[3]
                    
                gscores.append(float(score))
            
        if os.path.exists(self.path('DOCK_RMSD')):
            rmsds = []
            with open(self.path('DOCK_RMSD')) as fp:
                for line in fp:
                    line = line.strip().split(',')
                    if line[3] == '"RMSD"': continue
                    rmsds.append(float(line[3].strip('"')))
        else:
            rmsds = [None]*len(gscores)
        return gscores, emodels, rmsds

    def path(self, name, extras = {}):
        extras.update(self.params)
        return self.paths[name].format(**extras)

class LigandManager:
    def __init__(self, protein, params, paths):
        self.protein = protein
        self.params = {k:v for k, v in params.items()}
        self.params['protein'] = protein
        self.paths = paths

        # Initialize ligand info.
        self.pdb = self.read_pdb()
        self.chembl = self.read_chembl()
        self.prepped = self.prepped()

        # Set default structure.
        self.st = None
        if not os.path.exists(self.path('GRID_ROOT')): return
        self.grids = sorted([l for l in os.listdir(self.path('GRID_ROOT'))
                             if l[0] != '.'])
        if not self.grids: return
        
        self.st = self.grids[0]

        self.params['struct'] = self.st
        
        self.mcss = MCSSController(self)
        self.helpers = {}

    def read_pdb(self):
        pdb = {}
        if os.path.exists(self.path('PDB')):
            for _, row in pd.read_csv(self.path('PDB')).iterrows():
                name = '{}_lig'.format(row['PDB ID'])
                smiles = row['Ligand SMILES']
                affinity = 1
                pdb[name] = (smiles, affinity)
        return pdb

    def read_chembl(self):
        chembl = {}
        _chembl = []
        for csv in glob(self.path('CHEMBL')):
            _chembl += [pd.read_csv(csv)]
        if _chembl:
            for _, row in pd.concat(_chembl).iterrows():
                name = '{}_lig'.format(row['ligand_chembl_id'])
                smiles = row['canonical_smiles']
                affinity = row['standard_value']
                chembl[name] = (smiles, affinity)
        return chembl

    def prepped(self):
        if not os.path.exists(self.path('LIGANDS_ROOT')): return set([])
        return set([ligand for ligand in os.listdir(self.path('LIGANDS_ROOT'))
                    if os.path.exists(self.path('LIGANDS', {'ligand': ligand}))])

    def get_pdb(self, num=None):
        pdb = sorted(self.pdb.keys())
        if num is None:
            return pdb
        return pdb[:num]

    def get_chembl(self):
        return sorted(self.chembl.keys())

    def lookup(self, ligand):
        if ligand in self.pdb:
            return self.pdb[ligand]
        if ligand in self.chembl:
            return self.chembl[ligand]
        return None

    def get_structure(self, ligand):
        return SmilesStructure(self.lookup(ligand)[0]).get2dStructure()

    def get_xdocked_ligands(self, num):
        ligands = self.docked(self.get_pdb())[:num+1]
        self_docked = self.st+'_lig'
        if self_docked in ligands:
           ligands.remove(self_docked)
        elif ligands:
            ligands = ligands[:num]
        return ligands

    def docked(self, ligands, st=None):
        return [ligand for ligand in ligands
                if os.path.exists(self.path('DOCK_PV', {'ligand': ligand}))]

    def unique(self, ligands_st, query_st):
        for ligand_st in ligands_st:
            if ligand_st.isEquivalent(query_st):
                return False
        return True

    def diverse(self, ligands, query):
        for ligand in ligands:
            if self.mcss.get_mcss_size(ligand, query, compute=True) > 0.8:
                return False
        return True

    def _pick_helpers_best_affinity(self, query, num_chembl):
        sorted_helpers = sorted(self.get_chembl(), key=lambda x: self.chembl[x][-1])
        query_st = self.get_structure(query)
        helpers, helpers_st = [], []
        for ligand in sorted_helpers:
            ligand_st = self.get_structure(ligand)
            if self.unique(helpers_st+[query_st], ligand_st):
                helpers += [ligand]
                helpers_st += [ligand_st]
            if len(helpers) == num_chembl:
                break
        return helpers

    def _pick_helpers_best_affinity_diverse(self, query, num_chembl):
        sorted_helpers = sorted(self.get_chembl(), key=lambda x: self.chembl[x][-1])
        query_st = self.get_structure(query)
        helpers, helpers_st = [], []
        for ligand in sorted_helpers:
            ligand_st = self.get_structure(ligand)
            if self.unique(helpers_st+[query_st], ligand_st) and self.diverse(helpers, ligand):
                helpers += [ligand]
                helpers_st += [ligand_st]
            if len(helpers) == num_chembl:
                break
            print(len(helpers))
        return helpers

    def _pick_helpers_best_affinity_maybe_diverse(self, query, num_chembl):
        sorted_helpers = sorted(self.get_chembl(), key=lambda x: self.chembl[x][-1])
        query_st = self.get_structure(query)
        helpers, helpers_st = [], []
        for ligand in sorted_helpers:
            ligand_st = self.get_structure(ligand)
            if self.unique(helpers_st+[query_st], ligand_st) and self.diverse(helpers, ligand):
                helpers += [ligand]
                helpers_st += [ligand_st]
            if len(helpers) == num_chembl:
                break

        if len(helpers) < num_chembl:
            for ligand in sorted_helpers:
                if ligand in helpers: continue
                if self.unique(helpers_st+[query_st], ligand_st):
                    helpers += [ligand]
                    helpers_st += [ligand_st]
                if len(helpers) == num_chembl:
                    break
        return helpers

    def _pick_helpers_best_mcss(self, query, num_chembl):
        sorted_helpers = sorted(self.get_chembl(), key=lambda x: self.chembl[x][-1])
        sorted_helpers = self.mcss.sort_by_mcss(query, sorted_helpers)
        query_st = self.get_structure(query)
        helpers, helpers_st = [], []
        for ligand in sorted_helpers:
            ligand_st = self.get_structure(ligand)
            if self.unique(helpers_st+[query_st], ligand_st):
                helpers += [ligand]
                helpers_st += [ligand_st]
            if len(helpers) == num_chembl:
                break
        return helpers

    def pick_helpers(self, maxnum=20, num_chembl=25):
        os.system('mkdir -p {}'.format(self.path('HELPERS_ROOT')))
        options = {
            'best_affinity': self._pick_helpers_best_affinity,
            'best_mcss': self._pick_helpers_best_mcss,
            'best_affinity_diverse': self._pick_helpers_best_affinity_diverse,
            'best_affinity_maybe_diverse': self._pick_helpers_best_affinity_maybe_diverse
            }

        for option, function in options.items():
            path = self.path('HELPERS', {'helpers': option})
            if os.path.exists(path): continue
            print('picking chembl ligands', option)
            self.mcss.load_mcss()
            with open(path, 'w') as fp:
                for query in self.get_xdocked_ligands(maxnum):
                    print(query)
                    helpers = function(query, num_chembl)
                    helpers = ','.join(helpers)
                    fp.write('{}:{}\n'.format(query, helpers))

    def load_helpers(self):
        helpers = {}
        glob_pattern = self.path('HELPERS', {'helpers': '*'})
        re_pattern = self.path('HELPERS', {'helpers': '(.+)'})
        for fname in glob(glob_pattern):
            fname = re.match(re_pattern, fname).group(1)
            helpers[fname] = {}
            with open(self.path('HELPERS', {'helpers': fname})) as fp:
                for line in fp:
                    q, chembl = line.strip().split(':')
                    helpers[fname][q] = chembl.split(',')
        return helpers

    def get_helpers(self, query, fname, num=10, struct=None, randomize=False):
        if struct is None: struct = self.st
        helpers = self.load_helpers()
        print(list(helpers.keys()))
        helpers = helpers[fname][query]
        helpers = self.docked(helpers, struct)
        if randomize:
            # This is probably overkill, but I don't want random noise as
            # I'm optimizing parameters.
            np.random.seed(hash(self.protein) % (2**32 - 1))
            idx = np.random.permutation(len(helpers))
            helpers = [helpers[i] for i in idx]
        return helpers[:num]

    def path(self, name, extras = {}):
        extras.update(self.params)
        return self.paths[name].format(**extras)

class Protein:
    def __init__(self, protein, params, paths):
        """
        protein (str): name of protein, should be name of directory in root
            of data directory.
        data (str): path to the data directory.
        params (dict): various parameters.
        """
        self.protein = protein
        self.params = {k:v for k, v in params.items()}
        self.params['protein'] = protein
        self.paths = paths
        
        self.lm = LigandManager(protein, self.params, self.paths)

        if self.lm.st:
            self.docking = {self.lm.st: {}}
        else:
            self.docking = {}

    def load_docking(self, ligands, load_fp=False, load_crystal=False,
                     load_mcss=False, st=None):
        if st is None:
            st = self.lm.st

        if st not in self.docking:
            self.docking[st] = {}

        for ligand in ligands:
            params = {'struct': st}
            params.update(self.params)
            self.docking[st][ligand] = Ligand(ligand, params, self.paths)
            self.docking[st][ligand].load_poses(load_fp)

        if load_mcss:
            ligands = ligands+list(self.docking[st].keys())
            self.lm.mcss.load_rmsds(ligands, self.params['max_poses'])
