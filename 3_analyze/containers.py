"""
These classes keep track of
Poses, Ligands, Docking results for a given grid,
ligand names and preparation progress (LigandManager),
and ligands + docking results for a Protein.

Design decisions:
 - Parameterize by passing directory paths
   - Exception for struct, since needs to be specified
     in the docking file name
 - Open all files with absolute paths
 - Only load data if explicitly specified
 - Have default structure specified by LigandManager.
   always use this structure unless explicitly specified
   otherwise.
"""

import os
import sys
import numpy as np

from shared_paths import shared_paths
from parse_files import parse_glide_output, parse_fp_file
from parse_chembl import load_chembl_proc
from pick_helpers import load_helpers
from chembl_props import read_duplicates
from mcss_controller import MCSSController

class Pose:
    """
    Represents a single pose of a ligand.

    Parameterized by (protein, structure, ligand, rank).
    """
    def __init__(self, rmsd, gscore, emodel, fp):
        """
        rmsd (float): RMSD to crystallographic pose
        gscore (float): glide score
        emodel (float): emodel (from glide)
        fp ({(int, string): float}): dict mapping interactiontype, resname pairs to interaction score
        """
        self.rmsd = rmsd
        self.gscore = gscore
        self.emodel = emodel
        self.fp = fp

class Ligand:
    """
    Stores all poses of a given ligand docked to a given structure.

    This class defines file path convention for docking and ifp subdirectories.

    Parameterized by (protein, structure, ligand).
    """
    def __init__(self, dock_dir, fp_dir, ligand, struct):
        self.glide_path = '{}/{}-to-{}'.format(dock_dir, ligand, struct)
        self.fp_path = '{}/{}-to-{}.fp'.format(fp_dir, ligand, struct)
        self.crystal_fp_path = '{}/{}.fp'.format(fp_dir, ligand)

        self.poses = None
        self.crystal_pose = None

    def load_poses(self, load_fp):
        gscores, emodels, rmsds = parse_glide_output(self.glide_path)

        fps = {}
        if load_fp:  
            fps = parse_fp_file(self.fp_path)
            assert len(fps) >= min(100, len(rmsds)), \
                   ('missing for {}'.format(self.fp_path))
 
        assert len(gscores) == len(rmsds) == len(emodels), \
               ('Number of glide scores, rmsds, and emodels '
                'is not equal for {}'.format(self.glide_path))
        assert gscores == sorted(gscores), \
               'Glide scores are not correctly ordered for {}'.format(self.glide_path)

        self.poses = [Pose(rmsds[i], gscores[i], emodels[i], fps.get(i, {}), i)
                      for i in range(len(gscores))]

    def load_crystal_pose(self):
        try:
            fps = parse_fp_file(self.crystal_fp_path)
            self.crystal_pose = Pose(0, 0, 0, fps[0], 0)
        except IOError:
            pass

class Docking:
    """
    Stores docking results for a set of ligands to a single structure.

    This class defines path to docking and ifp directories.

    Parameterized by (protein, structure).
    """
    def __init__(self, root, struct):
        self.dock_dir = '{}/docking/{}'.format(root, shared_paths['docking'])
        self.ifp_dir = '{}/ifp/{}'.format(root, shared_paths['ifp'])
        self.struct = struct
        self.ligands = {}
        self.num_poses = {}

    def load(self, ligands, load_fp, load_crystal_fp):
        for ligand in ligands:
            pair = '{}-to-{}'.format(ligand, self.struct)
            self.ligands[l] = Ligand(ligand, self.dock_dir, self.ifp_dir, self.struct)
            self.ligands[l].load_poses(load_fp)
            if load_crystal_fp: self.ligands[l].load_crystal_pose()
            self.num_poses[l] = len(self.ligands[l].poses)

class LigandManager:
    """
    Manages docking results, fps, and MCSSs for all ligands associated
    with a given protein.

    This class defines path to protein directory.

    Parameterized by (protein).
    """
    def __init__(self, root, struct = None):
        """
        root (string): Path to protein's data directory.
        """
        self.root = root
        
        # Initialize ligand info.
        self.chembl_info = load_chembl_proc(self.root)
        self.u_ligs, self.dup_ligs = read_duplicates(self.root)
        self.all_ligs = self.prepped()
        self.pdb = self.unique(sorted([l for l in self.all_ligs if l[:6] != 'CHEMBL']))

        # Set default structure.
        self.st = None
        if not os.path.exists('{}/docking/grids'.format(self.root)): return
        self.grids = sorted([l for l in os.listdir('{}/docking/grids'.format(self.root)) if l[0] != '.'])
        if not self.grids: return
        self.st = self.grids[0] if struct is None else struct

        self.mcss = MCSSController(self)
        self.helpers = {}

    def docked(self, ligands, st=None):
        if st == None: st = self.st
        return [l for l in ligands
                if os.path.exists(Ligand(Docking(self.root, st).dock_dir,
                                         '', ligand, st).glide_path)]

    def prepped(self):
        ligdir = '{}/ligands/prepared_ligands'.format(self.root)
        if not os.path.exists(ligdir): return set([])
        return set([l for l in os.listdir(ligdir) if os.path.exists('{}/{}/{}_out.mae'.format(ligdir,l,l))])

    def chembl(self):
        def valid(l):
            filters = [
                lambda x,ci: ci[x].ki is not None and ci[x].ki <= 1000,
                lambda x,ci: ci[x].mw is not None and ci[x].mw <= 800,
                lambda x,ci: ci[x].macrocycle is not None and not ci[x].macrocycle
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

    def get_helpers(self, query, fname, num=10, struct=None):
        if struct is None: struct = self.st
        if fname not in self.helpers:
            self.helpers[fname] = load_helpers(self.root)[fname]
            for q in self.helpers[fname]:
                self.helpers[fname][q] = self.docked(self.helpers[fname][q], struct)
        return self.helpers[fname][query][:num]

class Protein:
    """
    Collection of ligands and docking results for a given protein.
    """
    def __init__(self, protein):
        self.root = "{}/{}".format(shared_paths['data'], protein)
        self.lm = LigandManager(self.root, None)
        self.docking = {}

    def load_docking(self, ligands, load_fp, load_crystal_fp, load_mcss, st = None):
        if st == None: st = self.lm.st

        if st not in self.docking:
            self.docking[st] = Docking(self.root, st)
        self.docking[st].load(ligands, load_fp, load_crystal_fp)

        if load_mcss:
            self.lm.mcss.load_rmsds(ligands, 100)
