"""
These classes keep track of
Poses, Ligands, Docking results for a given grid,
ligand names and preparation progress (LigandManager),
and ligands + docking results for a Protein.
"""

import os
import sys
import numpy as np

from ifp.fp_controller import parse_fp_file
from dock.parse_chembl import load_chembl_proc
from dock.chembl_props import read_duplicates
from mcss.mcss_controller import MCSSController

class Pose:
    def __init__(self, rmsd, gscore, emodel, fp):
        """
        rmsd (float): RMSD to crystallographic pose
        gscore (float): glide score
        emodel (float): emodel (from glide)
        fp ({(int, string): float}): dict mapping (interactiontype, resname)
            to interaction score
        """
        self.rmsd = rmsd
        self.gscore = gscore
        self.emodel = emodel
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

        self.poses = [Pose(rmsds[i], gscores[i], emodels[i], fps.get(i, {}))
                      for i in range(len(gscores))]

    def load_crystal_pose(self):
        fps = parse_fp_file(self.path('XTAL_IFP'))
        self.poses = [Pose(0, -1000, -1000, fps[0])]

    def parse_glide_output(self):
        if not os.path.exists(self.path('DOCK')):
            return [], [], []

        gscores, emodels, rmsds = [], [], []
        with open(self.path('DOCK_REPT')) as fp:
            for line in fp:
                line = line.strip().split()
                if len(line) <= 1: continue
                lig, _lig = self.params['ligand'], line[1]
                if not (_lig == lig or _lig == lig+'_out' or _lig == '1'):
                    continue
                
                if _lig == '1':
                    score = line[2]
                    emodel = line[12]
                else:
                    score = line[3]
                    emodel = line[13]
                    
                gscores.append(float(score))
                emodels.append(float(emodel))
            
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
        self.chembl_info = load_chembl_proc(self.path('ROOT'))
        self.u_ligs, self.dup_ligs = read_duplicates(self.path('ROOT'))
        self.all_ligs = self.prepped()
        self.pdb = self.unique(sorted([l for l in self.all_ligs
                                      if l[:6] != 'CHEMBL']))

        # Set default structure.
        self.st = None
        if not os.path.exists(self.path('GRID_ROOT')): return
        self.grids = sorted([l for l in os.listdir(self.path('GRID_ROOT'))
                             if l[0] != '.'])
        if not self.grids: return
        
        if self.params['pdb_order'] is 'First':
            self.st = self.grids[0] 
        elif self.params['pdb_order'] is 'Last':
            self.st = self.grids[-1]
        else:
            assert False
        self.params['struct'] = self.st
        
        self.mcss = MCSSController(self)
        self.helpers = {}

    def get_xdocked_ligands(self, num):
        ligands = self.docked(self.pdb)[:num+1]
        self_docked = self.st+'_lig'
        if self_docked in ligands:
           ligands.remove(self_docked)
        elif ligands:
            ligands.pop(-1)
        return ligands

    def docked(self, ligands, st=None):
        if st == None: st = self.st
        return [ligand for ligand in ligands
                if os.path.exists(self.path('DOCK_PV', {'ligand': ligand}))]

    def prepped(self):
        if not os.path.exists(self.path('PREPARED_ROOT')): return set([])
        return set([ligand for ligand in os.listdir(self.path('PREPARED_ROOT'))
                   if os.path.exists(self.path('PREPARED', {'ligand': ligand}))])

    def chembl(self):
        def valid(l):
            filters = [
                lambda x,ci: ci[x].ki is not None and ci[x].ki <= 1000,
                lambda x,ci: ci[x].mw is not None and ci[x].mw <= 800,
                lambda x,ci: ci[x].macrocycle is not None and not ci[x].macrocycle,
                lambda x,ci: ci[x].valid_stereo is not None and ci[x].valid_stereo
            ]
            return all(f(l, self.chembl_info) for f in filters)
        return [l for l in self.all_ligs if l in self.chembl_info and valid(l)]

    def unique(self, l_list):
        """
        Removes duplicates from l_list.
        
        If identical ligands are found, the one that appears first in l_list will be kept.
        """
        if not self.u_ligs:
            print('duplicates not loaded')
            return l_list
        unique_ligs = []
        exclude = set([])
        for l in l_list:
            if l in exclude: continue
            unique_ligs.append(l)
            if l in self.u_ligs: continue
            for l2, dups in self.dup_ligs.items():
                if l in dups:
                    exclude.update(dups)
                    break
            else:
                print('uh oh, ligand not found in unique or duplicates...', l)
        return unique_ligs

    def pick_helpers(self, maxnum=20, num_chembl=20):
        ki_sort = lambda c: self.chembl_info[c].ki

        # filters
        all_options = [
            'best_affinity.txt',
            'best_mcss.txt',
            'best_affinity_diverse.txt'
        ]

        root = self.path('HELPERS')
        os.system('mkdir -p {}'.format(root))
        for f in all_options:
            path = '{}/{}'.format(root, f)
            if os.path.exists(path): continue
            print('picking chembl ligands', f)
            chembl_ligs = sorted(self.chembl())
            with open(path, 'w') as fi:
                for q in self.get_xdocked_ligands(maxnum):
                    sorted_helpers = sorted(chembl_ligs, key=ki_sort)
                    print(q)
                    if f == 'best_mcss.txt':
                        self.mcss.load_mcss()
                        sorted_helpers = self.mcss.sort_by_mcss(q, sorted_helpers)

                    unique = self.unique([q]+sorted_helpers)
                    if f == 'best_affinity_diverse.txt':
                        self.mcss.load_mcss()
                        _unique = []
                        for helper in unique:
                            for _helper in _unique:
                                not_query = _helper != q
                                sim = self.mcss.get_mcss_size(helper, _helper,compute=True) > 0.8
                                if sim and not_query:
                                    break
                            else:
                                _unique += [helper]

                            print(len(_unique))
                            if len(_unique) == num_chembl+1:
                                break
                        unique = _unique
                    fi.write('{}:{}\n'.format(q, ','.join(unique[1:num_chembl+1])))

    def load_helpers(self):
        root = self.path('HELPERS')
        helpers = {}
        for fname in os.listdir(root):
            if fname[0] == '.' or fname.split('.')[-1] != 'txt': continue
            helpers[fname] = {}
            with open('{}/{}'.format(root, fname)) as fp:
                for line in fp:
                    q, chembl = line.strip().split(':')
                    helpers[fname][q] = chembl.split(',')
        return helpers

    def get_helpers(self, query, fname, num=10, struct=None, randomize=False):
        if struct is None: struct = self.st
        if fname not in self.helpers:
            self.helpers[fname] = self.load_helpers()[fname]
            for q in self.helpers[fname]:
                self.helpers[fname][q] = self.docked(self.helpers[fname][q], struct)
        if randomize:
            # This is probably overkill, but I don't want random noise as
            # I'm optimizing parameters.
            np.random.seed(hash(self.protein) % (2**32 - 1))
            idx = np.random.permutation(len(self.helpers[fname][query]))
            helpers = [self.helpers[fname][query][i] for i in idx]
        else:
            helpers = self.helpers[fname][query]
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
            if load_crystal:
                assert 'crystal' in ligand
                self.docking[st][ligand].load_crystal_pose()
            else:
                assert 'crystal' not in ligand
                self.docking[st][ligand].load_poses(load_fp)

        if load_mcss:
            ligands = ligands+list(self.docking[st].keys())
            self.lm.mcss.load_rmsds(ligands, self.params['max_poses'])
