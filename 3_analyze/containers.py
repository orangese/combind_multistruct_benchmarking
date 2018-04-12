import os
import sys

from parse_files import parse_glide_output, parse_fp_file, parse_mcss, parse_mcss_size
from plotting_tools import plot_docking

sys.path.append('../1_dock')
from parse_chembl import load_chembl_raw, load_chembl_proc
from ligand_features import load_features

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
        all_pairs = [(l1,l2) for l1 in ligands for l2 in ligands if (l1,l2) not in self.mcss]
        new_pairs = parse_mcss(self.mcss_dir, self.num_poses, all_pairs)
        for new in new_pairs:
            self.mcss[new] = new_pairs[new]

    def get_mcss_score(self, l1, p1, l2, p2):
        if (l1,l2) in self.mcss: return self.mcss[(l1,l2)][(p1,p2)]
        if (l2,l1) in self.mcss: return self.mcss[(l2,l1)][(p2,p1)]

        return None

    def load(self, load_fp, ligands):
        for l in ligands:
            if l in self.ligands: continue
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
    def __init__(self, uniprot, dock_st, data_dir, glide_dir, ifp_dir, mcss_dir):

        self.uniprot = uniprot
        self.dock_st = dock_st
        
        self.lm = LigandManager(uniprot, dock_st, data_dir, glide_dir, mcss_dir)

        self.dock_dir = '{}/{}/{}'.format(data_dir, uniprot, glide_dir)
        self.ifp_dir = '{}/{}/{}'.format(data_dir, uniprot, ifp_dir)
        self.mcss_dir = '{}/{}/{}/{}'.format(data_dir, uniprot, mcss_dir, dock_st)

        self.docking = Docking(self.dock_dir, self.ifp_dir, self.mcss_dir, self.dock_st)
        self.true = {}

    def load_true(self, ifp_dir, ligs):
        for l in ligs:
            self.true[l] = Ligand(l)
            self.true[l].load_poses(None, ifp_dir, l)

class LigandManager:
    def __init__(self, prot, struct, data_dir, gdir, mdir):
        self.root = '{}/{}'.format(data_dir, prot)
        self.struct = struct
        self.gdir = gdir
        self.mdir = mdir

        self.unique = set([l.split('.')[0] for l in os.listdir('{}/ligands/unique'.format(self.root))])
        self.chembl_info = {'{}_lig'.format(l):info for l,info in load_chembl_proc(self.root).items()}

        self.docked = []
        self.mcss_sizes = {}

    def init_docked(self):
        glide_dir = '{}/{}'.format(self.root, self.gdir)
        for fname in os.listdir(glide_dir):
            try: lig, st = fname.split('-to-')
            except: continue
            if st == self.struct and os.path.exists('{}/{}/{}_pv.maegz'.format(glide_dir, fname, fname)):
                self.docked.append(lig)
        self.docked.sort()

    def get_pdb(self):
        if len(self.docked) == 0:
            self.init_docked()
        return [l for l in self.docked if l[:6] != 'CHEMBL']

    def get_chembl(self, stereo=True, max_ki=float('inf')):
        if len(self.docked) == 0:
            self.init_docked()
        tr = [l for l in self.docked if l[:6] == 'CHEMBL' and (not stereo or self.chembl_info[l].valid_stereo)]
        tr = [l for l in tr if self.chembl_info[l].ki <= max_ki]
        return sorted(tr, key=lambda x: self.chembl_info[x].ki)

    def get_similar(self, query, num=10, stereo=True, max_ki=float('inf'), chembl=True):
        all_ligs = self.get_chembl(stereo, max_ki)
        if not chembl:
            all_ligs = [l for l in self.get_pdb() if l != query]

        if query not in self.mcss_sizes:
            mcss_path = '{}/{}/{}'.format(self.root, self.mdir, self.struct) 
            all_pairs = parse_mcss_size(mcss_path, set([query]), set(all_ligs))
            self.mcss_sizes[query] = {}
            for (l1,l2),sz in all_pairs.items():
                if query == l1: 
                    self.mcss_sizes[query][l2] = sz[2]
                if query == l2:
                    self.mcss_sizes[query][l1] = sz[2]

        return sorted(all_ligs, key=lambda x:-self.mcss_sizes[query].get(x,0))[:num]

class Dataset:
    def __init__(self, prots, structs, data_dir, glide_dir, ifp_dir, mcss_dir):
        self.data_dir = data_dir
        self.proteins = {}
        for u in prots:
            self.proteins[u] = Protein(u, structs[u], data_dir, glide_dir, ifp_dir, mcss_dir)
            

    def load(self, ligs={}, load_fp=True, load_crystal=False, load_mcss=True):
        for p, l_list in ligs.items():
            dock = self.proteins[p].docking
            dock.load(load_fp, l_list)

            if load_mcss:
                lig_feats = load_features('{}/{}'.format(self.data_dir, p))
                dock.load_mcss(l_list)
                for l, lig in dock.ligands.items():
                    lig.chrg = 'chrg' in lig_feats[l]
                    lig.ring = 'ring' in lig_feats[l]

            if load_crystal:
                self.proteins[p].load_true(load_fp)

    def assign_weights(self, w):
        for p_name, p in self.proteins.items():
            for l_name, l in p.docking.ligands.items():
                for pose in l.poses:
                    pose.weight_fp(w)
            for l_name, l in p.true.items():
                for pose in l.poses:
                    pose.weight_fp(w)

    def glide_perf(self, g_name, uniprots, ligands=None, show=True, n_list=[1,5,25,100], title=''):
        all_rmsds = [ [] for n in n_list ]
        for u in uniprots:
            p = self.proteins[u]
            rmsds = p.docking[g_name].glide_perf(n_list=n_list, ligands=ligands)
            for i, n in enumerate(n_list):
                all_rmsds[i].extend(rmsds[i])
        if show and len(rmsds[0]) > 0:
            from plotting_tools import plot_docking
            plot_docking(all_rmsds, ['top-{}'.format(n) for n in n_list], title)
        return all_rmsds

