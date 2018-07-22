import os
import sys
import numpy as np

from parse_files import parse_glide_output, parse_fp_file

from mcss_utils import MCSS

sys.path.append('../1_dock')
sys.path.append('/'.join(os.path.dirname(os.path.realpath(sys.argv[0])).split('/')[:-1]) + '/1_dock')

from parse_chembl import load_chembl_proc
from pick_helpers import load_helpers
from chembl_props import read_duplicates

class Pose:
    """
    Represents a single pose of a ligand.
    """
    def __init__(self, rmsd, physics_score, fp, rank):
        """
        rmsd (float): RMSD to crystallographic pose
        physics_score (float): score from physics-based docking alg
        fp ({(int, string): float}): dict mapping interactiontype, resname pairs to interaction score
        rank (int): ranking by physics score, and thereby index in pose file
        """
        self.rmsd = rmsd
        self.gscore = physics_score
        self.fp_raw = fp
        self.rank = rank

        self.fp = fp

    def weight_fp(self, w):
        self.fp = {}
        for (i, r), sc in self.fp_raw.items():
            if i in w and w[i] != 0:
                self.fp[(i,r)] = sc*w[i]

class Ligand:
    """
    Stores all poses of a given ligand.
    """
    def __init__(self, lig_id):

        self.lig_id = lig_id
        self.poses = []

    def load_poses(self, dock_dir, fp_dir, struct, load_fp):
        if struct == self.lig_id:
            fp_path = '{}/{}.fp'.format(fp_dir, self.lig_id)
            fps = parse_fp_file(fp_path)
            try:
                self.poses = [Pose(0, 0, fps[0], 'T')]
            except: print(fp_path)
            return

        f_path = '{}/{}-to-{}'.format(dock_dir, self.lig_id, struct)
        gscores, rmsds = parse_glide_output(f_path)

        fp_path = '{}/{}-to-{}.fp'.format(fp_dir, self.lig_id, struct)
        fps = {}
        if load_fp:
            fps = parse_fp_file(fp_path)
            if len(fps) < min(100, len(rmsds)):
                print('missing fp?', fp_path)
 
        if not len(gscores) == len(rmsds):
            print('{} {} {}'.format(dock_dir, self.lig_id, struct))
            print(f_path, fp_path)
            #os.system('rm -rf {}'.format(f_path))
            #os.system('rm -rf {}'.format(fp_path))
            self.poses = []
            return            

        self.poses = [Pose(rmsds[i], gscores[i], fps.get(i, {}), i) for i in range(len(gscores))]
        
        for i, p in enumerate(self.poses):
            if i == 0: continue
            assert p.gscore >= self.poses[i-1].gscore, '{} {} {}'.format(dock_dir, self.lig_id, struct)
                
        if len(self.poses) == 0: print(len(rmsds), len(gscores), len(fps))

    def top_n(self, n):
        if len(self.poses) == 0:
            return Pose(100, 0, {}, 300)
        return self.poses[np.argmin([self.poses[i].rmsd for i in range(min(n, len(self.poses)))])]

class Docking:
    """
    Stores docking results for a set of ligands.
    """
    def __init__(self, sp, prot, struct):
        self.sp = sp
        self.dock_dir = '{}/{}/docking/{}'.format(sp['data'], prot, sp['docking'])
        self.ifp_dir = '{}/{}/ifp/{}'.format(sp['data'], prot, sp['ifp'])

        self.struct = struct
        self.ligands = {}
        self.num_poses = {}

    def load(self, ligands, load_fp):
        for l in ligands:
            if l in self.ligands: 
                if not (load_fp and self.ligands[l].poses[0].fp == {}):
                    continue
            pair = '{}-to-{}'.format(l, self.struct)
            if not os.path.exists('{}/{}/{}_pv.maegz'.format(self.dock_dir, pair, pair)): continue
            self.ligands[l] = Ligand(l)
            self.ligands[l].load_poses(self.dock_dir, self.ifp_dir, self.struct, load_fp)
            self.num_poses[l] = len(self.ligands[l].poses)

    def glide_perf(self, n_list=[1,5,25,100], ligands=None):
        rmsds = [ [] for n in n_list ]
        if ligands is None:
            ligands = self.ligands.keys()
        for i,n in enumerate(n_list):
            for l in ligands:
                if l in self.ligands:
                    rmsds[i].append(self.ligands[l].top_n(n).rmsd)
                else:
                    rmsds[i].append(None)
        return rmsds

class Protein:
    """
    Stores ligands and docking results for a given protein.
    """
    def __init__(self, shared_paths, prot, struct):
        self.prot = prot
        self.lm = LigandManager(shared_paths, prot, struct)
        self.sp = shared_paths

        # if a struct is provided (above), lm.st will use it
        # otherwise lm.st will provide a default
        self.docking = { self.lm.st : Docking( shared_paths, self.prot, self.lm.st) }
        self.lm.mcss.num_poses = self.docking[self.lm.st].num_poses
        
        self.true = {}

    def load(self, l_list, st, load_fp, load_crystal, load_mcss):
        if st == None: st = self.lm.st

        if st not in self.docking:
            self.docking[st] = Docking(self.sp, self.prot, st)
        self.docking[st].load(l_list, load_fp)

        if load_mcss:
            self.lm.mcss.load_mcss(l_list, l_list, rmsd=True)        

        if load_crystal:
            for l in l_list:
                if l in self.true: continue
                self.true[l] = Ligand(l)
                self.true[l].load_poses(None, ifp_dir, l)

class LigandManager:
    """
    Manages docking results, fps, and MCSSs for all ligands associated
    with a given target.
    """
    def __init__(self, shared_paths, prot, struct=None):
        """
        shared_paths ({string: string}): dict at least containing the key
                                         'data' providing the path to where
                                         data is stored.
        prot (string): the name of the desired target.
        """
        # Get path to data
        self.root = '{}/{}'.format(shared_paths['data'], prot)
        self.prot = prot
        self.sp = shared_paths

        self.all_st = {} # What is this for? Vestigial?
        self.st = None

        # Get ligand info
        self.chembl_info = load_chembl_proc(self.root)
        self.u_ligs, self.dup_ligs = read_duplicates(self.root)
        self.all_ligs = self.prepped()
        self.pdb = self.unique(sorted([l for l in self.all_ligs if l[:6] != 'CHEMBL']))

        # As written, this will always get the first structure.
        if not (os.path.exists('{}/docking/grids'.format(self.root))
                and os.listdir('{}/docking/grids'.format(self.root))): return

        self.grids = sorted([l for l in os.listdir('{}/docking/grids'.format(self.root)) if l[0] != '.'])
        self.first_st = self.grids[0]
        self.st = struct
        if struct is None:
            self.st = self.all_st.get(prot, self.first_st)

        # Load MCSS
        self.mcss = MCSS(self.sp, self.st, self.root)
        self.helpers = {}

    def prepped(self):
        ligdir = '{}/ligands/prepared_ligands'.format(self.root)
        if not os.path.exists(ligdir): return set([])
        return set([l for l in os.listdir(ligdir) if os.path.exists('{}/{}/{}_out.mae'.format(ligdir,l,l))])

    def chembl(self, filters=[], sort_key=lambda x:0, unique=False):
        default = [
            lambda x,ci: ci[x].ki <= 1000,
            lambda x,ci: ci[x].mw <= 1000
        ]
        c = sorted([l for l in self.all_ligs if l in self.chembl_info 
                    and False not in [f(l,self.chembl_info) for f in default+filters]], key=sort_key)
        if unique: 
            return self.unique(c)
        return c

    def unique(self, l_list):
        # removes duplicates from l_list
        # if identical ligands are found, the one that
        # appears first in l_list will be kept
        if len(self.u_ligs) == 0: 
            'duplicates not loaded'
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

    def docked(self, l_list, st=None):
        if st == None: st = self.st
        gpath = '{}/docking/{}/{}-to-{}/{}-to-{}_pv.maegz'
        return [l for l in l_list if os.path.exists(gpath.format(self.root, self.sp['docking'], l,st,l,st))]

    def get_similar(self, query, fname, num=10, struct=None):
        if struct is None: struct = self.st
        if fname not in self.helpers:
            self.helpers[fname] = load_helpers(self.root)[fname]
            for q in self.helpers[fname]:
                self.helpers[fname][q] = self.docked(self.helpers[fname][q], struct)

        return self.helpers[fname][query][:num]

class Dataset:
    """
    A set of proteins.
    """
    def __init__(self, shared_paths, prots, structs={}):
        self.sp = shared_paths
        self.proteins = { p : Protein( shared_paths, p, structs.get(p,None) ) for p in prots }
            
    def load(self, ligs={}, structs={}, load_fp=True, load_crystal=False, load_mcss=True):
        for p, l_list in ligs.items():
            for st in structs.get(p, [None]):
                dock = self.proteins[p].load(l_list, st, load_fp, load_crystal, load_mcss)

    def assign_weights(self, w):
        for p_name, p in self.proteins.items():
            for st_name, st in p.docking.items():
                for l_name, l in st.ligands.items():
                    for pose in l.poses:
                        pose.weight_fp(w)
            for l_name, l in p.true.items():
                for pose in l.poses:
                    pose.weight_fp(w)
