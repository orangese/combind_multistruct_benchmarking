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
from ifp.fp_controller import parse_fp_file
from dock.parse_chembl import load_chembl_proc
from dock.pick_helpers import load_helpers
from dock.chembl_props import read_duplicates
from mcss.mcss_controller import MCSSController

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
    def __init__(self, ligand, dock_dir, fp_dir, struct):
        self.ligand = ligand
        self.glide_path = '{}/{}-to-{}'.format(dock_dir, ligand, struct)
        self.fp_path = '{}/{}-to-{}-{}.fp'.format(fp_dir, ligand, struct,
                                                  shared_paths['docking'])
        self.crystal_fp_path = '{}/{}.fp'.format(fp_dir,
                                                 ligand.replace('_crystal_lig', '_struct'))

        self.poses = None

    def load_poses(self, load_fp):
        gscores, emodels, rmsds = self.parse_glide_output()
        print(rmsds)

        fps = {}
        if load_fp:  
            fps = parse_fp_file(self.fp_path)
            assert len(fps) >= min(shared_paths['stats']['max_poses'], len(rmsds)), \
                   ('missing for {}'.format(self.fp_path))
 
        assert len(gscores) == len(rmsds) == len(emodels), \
               ('Number of glide scores, rmsds, and emodels '
                'is not equal for {}'.format(self.glide_path))
        assert gscores == sorted(gscores), \
               'Glide scores are not correctly ordered for {}'.format(self.glide_path)

        self.poses = [Pose(rmsds[i], gscores[i], emodels[i], fps.get(i, {}))
                      for i in range(len(gscores))]

    def load_crystal_pose(self):
        try:
            fps = parse_fp_file(self.crystal_fp_path)
            self.poses = [Pose(0, 0, 0, fps[0])]
        except IOError:
            pass

    def parse_glide_output(self):
        if not os.path.exists(self.glide_path):
            print("Glide path not found!!!")
            return [], [], []
        print("Glide path found!!!")
        pair = self.glide_path.split('/')[-1]
        if os.path.exists('{}/{}_pv.maegz'.format(self.glide_path, pair)):
            print("Rept path found!!!")
            return self.parse_rept_file(pair)
        else:
            print("Rept path not found!!!")
            print('not finished', self.glide_path)
            return [], [], []

    def parse_rept_file(self, pair):
        lig, prot = pair.split('-to-')
        rept_file = '{}/{}.rept'.format(self.glide_path, pair)
        rmsd_file = '{}/rmsd.csv'.format(self.glide_path, pair)
        print("RMSD filename is {}".format(rmsd_file))
        
        gscores, emodels, rmsds = [], [], []
        with open(rept_file) as fp:
            for line in fp:
                line = line.strip().split()
                if len(line) <= 1 or (line[1] != lig and line[1] != lig+'_out' and line[1] != '1'): continue
                rank, lig_name, lig_index, score = line[:4]
                emodel = line[13]
                if line[1] == '1':
                    rank, lig_index, score = line[:3]
                    emodel = line[12]
                gscores.append(float(score))
                emodels.append(float(emodel))
            
        if not os.path.exists(rmsd_file):
            return gscores, emodels, [None]*len(gscores)

        with open(rmsd_file) as fp:
            for line in fp:
                line = line.strip().split(',')
                if line[3] == '"RMSD"': continue
                rmsds.append(float(line[3][1:-1]))

        return gscores, emodels, rmsds

class Docking:
    """
    Stores docking results for a set of ligands to a single structure.

    This class defines path to docking and ifp directories.

    Parameterized by (protein, structure).
    """
    def __init__(self, root, struct):
        self.dock_dir = '{}/docking/{}'.format(root, shared_paths['docking'])
        self.ifp_dir = '{}/ifp/{}'.format(root, shared_paths['ifp']['version'])
        self.struct = struct
        self.ligands = {}
        self.num_poses = {}

    def load(self, ligands, load_fp, load_crystal):
        for ligand in ligands:
            pair = '{}-to-{}'.format(ligand.replace('_crystal', ''), self.struct)
            self.ligands[ligand] = Ligand(ligand, self.dock_dir, self.ifp_dir, self.struct)
            if load_crystal:
                self.ligands[ligand].load_crystal_pose()
            else:
                print("Loading poses for {}".format(ligand))
                self.ligands[ligand].load_poses(load_fp)
            self.num_poses[ligand] = len(self.ligands[ligand].poses)

class LigandManager:
    """
    Manages docking results, fps, and MCSSs for all ligands associated
    with a given protein.

    This class defines path to protein directory.

    Parameterized by (protein).
    """
    def __init__(self, protein, root, struct = 'First'):
        """
        root (string): Path to protein's data directory.
        """
        self.protein = protein
        self.root = root
        self.read_root = "{}/{}".format(shared_paths['read_data'], protein)
        self.write_root = "{}/{}".format(shared_paths['write_data'], protein)
        
        # Initialize ligand info.

        # Get a dict from ligand_id -> CHEMBL object representing ligand. This dict contains all chembl
        # ligands (i.e. ligands in chembl_info.tx) Note: this includes DUD-E ligands
        self.chembl_info = load_chembl_proc(self.write_root)

        # Get a set of unique ligands and a dict of duplicates, where the key is the first duplicate
        # (the numbering is arbitrary, "first" doesn't really mean anything) and the value is a set of
        # ligands that are duplicates of each other
        self.u_ligs, self.dup_ligs = read_duplicates(self.root)

        # Get a list of all ligands that are in ligands/prepared_ligands/
        self.all_ligs = self.prepped()

        # TODO: this doesn't seem like it'll work anymore
        ''''
        self.pdb = self.unique(sorted([l for l in self.all_ligs
                                      if l[:6] != 'CHEMBL'],
                                      reverse = struct == 'Last'))
                                      '''

        # Get the list composed of all unique PDB ligands_ids (str) in ligands/prepared_ligands/
        # I chose to identify a PDB ligand as a ligand s.t. there are exactly 4 characters before
        # the first '_'

        # print([ligand for ligand in self.all_ligs if len(ligand.split('_')[0]) == 4])
        self.pdb = self.unique(sorted([l for l in self.all_ligs
                                      if len(l.split('_')[0]) == 4],
                                      reverse = struct == 'Last'))
        # print(self.pdb)

        # Set the value of self.st as the default PDB structure. If this is the first struct, choose
        # the alphanumerically first PDB structure. If this is the last, choose the last PDB structure
        self.st = None
        if not os.path.exists('{}/docking/grids'.format(self.read_root)): return
        self.grids = sorted([l for l in os.listdir('{}/docking/grids'.format(self.read_root)) if l[0] != '.'])
        if not self.grids: return
        if struct is 'First':
            self.st = self.grids[0] 
        elif struct is 'Last':
            self.st = self.grids[-1]
        else:
            assert False

        self.mcss = MCSSController(self)
        self.helpers = {}

    def get_xdocked_ligands(self, num):
        ''' Get the first num PDB ligands, leaving out the ligand bound to the default PDB structure
        
        Returns:
        * list of num ligand_ids
        '''
        def my_docked(ligands, st=None):
            if st == None: st = self.st
            return [ligand for ligand in ligands
                    if os.path.exists("{}/{}-to-{}_pv.maegz".format(
                                    Ligand(ligand, Docking(self.read_root, st).dock_dir,
                                            '', st).glide_path,
                                    ligand, st
                                  ))]

        ligands = my_docked(self.pdb)[:num+1]
        if len(ligands) == 0: return []
        self_docked = self.st+'_lig'
        if self_docked in ligands:
           ligands.remove(self_docked)
        else:
            ligands.pop(-1)
        # print("In xdocked", ligands)
        return ligands

    def docked(self, ligands, st=None):
        if st == None: st = self.st
        return [ligand for ligand in ligands
                if os.path.exists("{}/{}-to-{}_pv.maegz".format(
                                  Ligand(ligand, Docking(self.root, st).dock_dir,
                                         '', st).glide_path,
                                  ligand, st
                                  ))]

    def prepped(self):
        ''' Get a set of ligand_id's of ligands that have been run through epik and prepwizard (i.e.
        ligands s.t. there exist a file data_root/ligands/prepared_ligands/ligand_id/ligand_id_out.mae)

        Returns:
        * set of ligand_id's
        '''
        ligdir = '{}/ligands/prepared_ligands'.format(self.write_root)
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
    def __init__(self, protein, struct = 'First'):
        self.root = "{}/{}".format(shared_paths['write_data'], protein)

        # Unless otherwise specified, instantiate a LigandManager object with the struct 'First'
        self.lm = LigandManager(protein, self.root, struct)

        # Useful to be able to reference this before loading data.
        self.docking = {self.lm.st: Docking(self.root, self.lm.st)}

    def load_docking(self, ligands, load_fp=False, load_crystal = False,
                     load_mcss = False, st = None):
        if load_crystal:
            assert all('crystal' in ligand for ligand in ligands)
        else:
            assert not any('crystal' in ligand for ligand in ligands)

        if st == None: st = self.lm.st

        if st not in self.docking:
            self.docking[st] = Docking(self.root, st)
        self.docking[st].load(ligands, load_fp, load_crystal)

        if load_mcss:
            self.lm.mcss.load_rmsds(ligands+list(self.docking[st].ligands.keys()),
                                    shared_paths['stats']['max_poses'])
