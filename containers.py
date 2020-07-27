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
from mcss.mcss_controller import MCSSController
from shape.shape_controller import ShapeController

class Pose:
    def __init__(self, rank, rmsd, gscore, fp):
        """
        rank (int): rank in original docking output.
        rmsd (float): RMSD to crystallographic pose
        gscore (float): glide score
        fp ({(int, string): float}): dict mapping (interactiontype, resname)
            to interaction score
        """
        self.rank = rank
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

    def load_native_poses(self, load_fp, thresh):
        self.load_poses(load_fp)
        self.poses = [pose for pose in self.poses if pose.rmsd <= thresh]

    def load_poses(self, load_fp):
        gscores, rmsds = self.parse_glide_output()

        fps = {}
        if load_fp:  
            fps = self.parse_ifp_file()

        if len(rmsds) < len(gscores):
            assert False, 'Not all RMSDs calculated for {}'.format(self.ligand)

        self.poses = [Pose(i, rmsds[i], gscores[i], fps.get(i, {}))
                      for i in range(len(gscores))]

    def parse_glide_output(self):
        if not os.path.exists(self.path('DOCK')):
            return [], []

        gscores = []
        if os.path.exists(self.path('DOCK_REPT')):
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

        rmsds = []
        if os.path.exists(self.path('DOCK_RMSD')):
            with open(self.path('DOCK_RMSD')) as fp:
                for line in fp:
                    line = line.strip().split(',')
                    if line[3] == '"RMSD"': continue
                    rmsds.append(float(line[3].strip('"')))
        
        if not rmsds and gscores:
            rmsds = [None]*len(gscores)
        elif rmsds and not gscores:
            gscores = [None]*len(rmsds)
        return gscores, rmsds

    def parse_ifp_file(self):
        ifps = {}
        if os.path.exists(self.path('IFP')):
            df = pd.read_csv(self.path('IFP'))
            for _, row in df.iterrows():
                pose = row['pose']
                res = row['protein_res']
                i = row['label']
                
                if i == 'saltbridge':
                    i = 1
                elif i == 'hbond_donor':
                    i = 2
                elif i == 'hbond_acceptor':
                    i = 3
                elif i == 'contact':
                    i = 11
                else:
                    assert False, i

                if pose not in ifps:
                    ifps[pose] = {}
                ifps[pose][(i,res)] = row['score']

        else:
            print(self.path('IFP'), 'fp not found')
        if len(ifps) == 0:
            print('check', self.path('IFP'))
            return {}
        return ifps

    def path(self, name, extras = {}):
        extras.update(self.params)
        return self.paths[name].format(**extras)

class LigandManager:
    def __init__(self, protein, params, paths, struct=None):
        self.protein = protein
        self.params = {k:v for k, v in params.items()}
        self.params['protein'] = protein
        self.paths = paths

        # Initialize ligand info.
        self.pdb = self.read_pdb()
        self.prepped = self.prepped()

        # Find available docking grids
        grids = []
        if os.path.exists(self.path('GRID_ROOT')):
            grids = sorted([l for l in os.listdir(self.path('GRID_ROOT'))
                            if l[0] != '.'])
        
        # Set default docking grid.
        if struct is not None:
            self.st = struct
        elif grids:
            self.st = grids[0]
        else:
            self.st = None

        self.params['struct'] = self.st

        if self.st in grids:
            self.mcss = MCSSController(self)
            self.shape = ShapeController(self)

    def add_ligands(self, fname):
        self.pdb.update(self.read_pdb(fname))

    def read_pdb(self, fname=None):
        if fname == None:
            fname = self.path('PDB')

        pdb = {}
        if os.path.exists(fname):
            for _, row in pd.read_csv(fname).iterrows():
                name = '{}_lig'.format(row['ID'])
                assert 'SMILES' in row
                pdb[name] = row
        return pdb

    def prepped(self):
        if not os.path.exists(self.path('LIGANDS_ROOT')): return set([])
        return set([ligand for ligand in os.listdir(self.path('LIGANDS_ROOT'))
                    if os.path.exists(self.path('LIGANDS', {'ligand': ligand}))])

    def get_pdb(self, num=None):
        pdb = sorted(self.pdb.keys())
        if num is None:
            return pdb
        return pdb[:num]

    def lookup(self, ligand):
        if ligand in self.pdb:
            return self.pdb[ligand]
        return None

    def get_xdocked_ligands(self, num=1000):
        ligands = self.docked(self.get_pdb())[:num+1]
        self_docked = self.st+'_lig'
        if self_docked in ligands:
           ligands.remove(self_docked)
        elif ligands:
            ligands = ligands[:num]
        return ligands

    def docked(self, ligands):
        return [ligand for ligand in ligands
                if os.path.exists(self.path('DOCK_PV', {'ligand': ligand}))]

    def path(self, name, extras = {}):
        extras.update(self.params)
        return self.paths[name].format(**extras)

class Protein:
    def __init__(self, protein, params, paths, struct=None):
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
        self.lm = LigandManager(protein, self.params, self.paths, struct)
        
        self.docking = {}

    def load_docking(self, ligands, load_fp=False, load_mcss=False, load_shape=False):
        for ligand in ligands:
            params = {'struct': self.lm.st}
            params.update(self.params)
            self.docking[ligand] = Ligand(ligand, params, self.paths)
            self.docking[ligand].load_poses(load_fp)

        if load_mcss:
            ligands = ligands+list(self.docking.keys())
            self.lm.mcss.load_rmsds(ligands, self.params['max_poses'])

        if load_shape:
            ligands = ligands+list(self.docking.keys())
            self.lm.shape.load(ligands)
