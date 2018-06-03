import os
import sys

from parse_files import parse_glide_output, parse_fp_file, parse_mcss, parse_mcss_size

sys.path.append('../1_dock')
from parse_chembl import load_chembl_raw, load_chembl_proc
from pick_helpers import load_helpers
from chembl_props import read_duplicates

import numpy as np

class Pose:
    def __init__(self, rmsd, physics_score, fp, rank):
        self.rmsd = rmsd
        self.gscore = physics_score
        self.fp_raw = fp
        self.fp = fp
        self.rank = rank

    def weight_fp(self, w):
        self.fp = {}
        for (i, r), sc in self.fp_raw.items():
            if i in w and w[i] != 0:
                self.fp[(i,r)] = sc*w[i]

class Ligand:
    def __init__(self, lig_id):
        self.lig_id = lig_id
        self.poses = []

    def load_poses(self, dock_dir, fp_dir, struct, load_fp):
        if struct == self.lig_id:
            fp_path = '{}/{}.fp'.format(fp_dir, self.lig_id)
            fps = parse_fp_file(fp_path)
            try:
                self.poses = [Pose(0, 0, fps[0], 'T')]
            except: print fp_path
            return

        f_path = '{}/{}-to-{}'.format(dock_dir, self.lig_id, struct)
        gscores, rmsds = parse_glide_output(f_path)

        fp_path = '{}/{}-to-{}.fp'.format(fp_dir, self.lig_id, struct)
        fps = {}
        if load_fp:
            fps = parse_fp_file(fp_path)
            if len(fps) < min(100, len(rmsds)):
                print 'missing fp?', fp_path       
 
        if not len(gscores) == len(rmsds):
            print '{} {} {}'.format(dock_dir, self.lig_id, struct)
            print f_path, fp_path
            os.system('rm -rf {}'.format(f_path))
            os.system('rm -rf {}'.format(fp_path))
            self.poses = []
            return            

        self.poses = [Pose(rmsds[i], gscores[i], fps.get(i, {}), i) for i in range(len(gscores))]
        
        for i, p in enumerate(self.poses):
            if i == 0: continue
            assert p.gscore >= self.poses[i-1].gscore, '{} {} {}'.format(dock_dir, self.lig_id, struct)
                
        if len(self.poses) == 0: print len(rmsds), len(gscores), len(fps)

    def top_n(self, n):
        if len(self.poses) == 0:
            return Pose(100, 0, {}, 300)
        return self.poses[np.argmin([self.poses[i].rmsd for i in range(min(n, len(self.poses)))])]

class Docking:
    def __init__(self, dock_dir, ifp_dir, mcss_dir, struct):
        self.dock_dir = dock_dir
        self.ifp_dir = ifp_dir
        self.mcss_dir = mcss_dir

        self.struct = struct
        self.ligands = {}
        self.num_poses = {}
        self.mcss = {}

    def load_mcss(self, ligands):
        all_pairs = [(l1,l2,self.struct) for l1 in ligands for l2 in ligands 
                        if (l1,l2) not in self.mcss and (l2,l1) not in self.mcss]
        new_pairs = parse_mcss(self.mcss_dir, self.num_poses, all_pairs)
        for new in new_pairs:
            self.mcss[new] = new_pairs[new]

    def get_mcss_score(self, l1, p1, l2, p2):
        if (l1,l2) in self.mcss: return self.mcss[(l1,l2)][(p1,p2)]
        if (l2,l1) in self.mcss: return self.mcss[(l2,l1)][(p2,p1)]

        return None

    def load(self, ligands, load_fp, load_mcss):
        for l in ligands:
            if l in self.ligands: continue
            pair = '{}-to-{}'.format(l, self.struct)
            if not os.path.exists('{}/{}/{}_pv.maegz'.format(self.dock_dir, pair, pair)): continue
            self.ligands[l] = Ligand(l)
            self.ligands[l].load_poses(self.dock_dir, self.ifp_dir, self.struct, load_fp)
            self.num_poses[l] = len(self.ligands[l].poses)
        if load_mcss:
            self.load_mcss(ligands)

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
    def __init__(self, uniprot, data_dir, glide_dir, ifp_dir, mcss_dir):
        self.uniprot = uniprot

        self.dock_dir = '{}/{}/{}'.format(data_dir, uniprot, glide_dir)
        self.ifp_dir = '{}/{}/{}'.format(data_dir, uniprot, ifp_dir)
        self.mcss_dir = '{}/{}/{}'.format(data_dir, uniprot, mcss_dir)

        self.lm = LigandManager(uniprot, data_dir, glide_dir, mcss_dir)

        self.true = {}
        self.docking = {}

    def load(self, l_list, st, load_fp, load_crystal, load_mcss):
        if st == None: st = self.lm.st

        if st not in self.docking:
            self.docking[st] = Docking(self.dock_dir, self.ifp_dir, self.mcss_dir, st)
        self.docking[st].load(l_list, load_fp, load_mcss)
        
        if load_crystal:
            for l in l_list:
                if l in self.true: continue
                self.true[l] = Ligand(l)
                self.true[l].load_poses(None, ifp_dir, l)

class LigandManager:
    def __init__(self, prot, data_dir, gdir, mdir):
        self.root = '{}/{}'.format(data_dir, prot)
        self.prot = prot
        self.gdir = gdir
        self.mdir = mdir

        #grids = {'D2R':'6CM4','AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
        #         'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23','TRPV1':'3J5Q}

        self.all_st = {}

        self.chembl_info = load_chembl_proc(self.root)
        self.u_ligs, self.dup_ligs = read_duplicates(self.root)

        self.all_ligs = self.prepped()
        self.pdb = self.unique(sorted([l for l in self.all_ligs if l[:6] != 'CHEMBL']))
        self.grids = sorted(os.listdir('{}/docking/grids'.format(self.root)))

        self.first_st = self.grids[0]
        self.st = self.all_st.get(prot, self.first_st)

        self.mcss_sizes = {}
        self.helpers = {}

    def prepped(self):
        ligdir = '{}/ligands/prepared_ligands'.format(self.root)
        if not os.path.exists(ligdir): return set([])
        return set([l for l in os.listdir(ligdir) if os.path.exists('{}/{}/{}.mae'.format(ligdir,l,l))])

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
                print 'uh oh, ligand not found in unique or duplicates...', l
        return unique_ligs

    def docked(self, l_list, st=None):
        if st == None: st = self.st
        gpath = '{}/{}/{}-to-{}/{}-to-{}_pv.maegz'
        return [l for l in l_list if os.path.exists(gpath.format(self.root, self.gdir, l,st,l,st))]

    def get_similar(self, query, fname, num=10, mcss_sort=False, struct=None):
        if struct is None: struct = self.st
        if fname not in self.helpers:
            self.helpers[fname] = load_helpers(self.root)[fname]
            for q in self.helpers[fname]:
                self.helpers[fname][q] = [l for l in self.helpers[fname][q] 
                    if l in set(self.get_docked(struct))]

        # mcss has two steps:
        # 1. we use the canvasMCS tool to find the MCSS between two ligands (output: smarts)
        # 2. we use the schrodinger python api to extract this mcss for each pose pair and compute the RMSD
        #    - to compute the rmsd, the two MCSS's must be "equivalent" according to the python api
        #    - some small fraction of the outputs from step 1 are not "equivalent" for unclear reasons
        #    - fname is sorted based on the output of step 1
        #    - here we re-sort based on the output of step 2 (identical to step 1 with failed pairs excluded)
        if mcss_sort:
            if len(self.chembl_info) == 0:
                self.chembl_info = load_chembl_proc(self.root, load_mcss=True)
            if fname not in self.mcss_sizes:
                self.mcss_sizes[fname] = {}
            if query not in self.mcss_sizes[fname]:
                mcss_path = '{}/{}'.format(self.root, self.mdir) 
                all_pairs = parse_mcss_size(mcss_path, set([query]), set(self.helpers[fname][query]), struct)
                self.mcss_sizes[fname][query] = {}
                for (l1,l2),sz in all_pairs.items():
                    if query == l1: 
                        self.mcss_sizes[fname][query][l2] = sz[2]
                    if query == l2:
                        self.mcss_sizes[fname][query][l1] = sz[2]

                ligs = self.helpers[fname][query]

                self.helpers[fname][query].sort(key=lambda x:-self.mcss_sizes[fname][query].get(x,0))

        return self.helpers[fname][query][:num]

class Dataset:
    def __init__(self, prots, data_dir, glide_dir, ifp_dir, mcss_dir):
        self.data_dir = data_dir
        self.proteins = {}
        for u in prots:
            self.proteins[u] = Protein(u, data_dir, glide_dir, ifp_dir, mcss_dir)
            
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


