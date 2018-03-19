import os
import sys

from parse_files import parse_glide_output, parse_fp_file, parse_mcss
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
            #print 'warning! {} did not dock'.format(self.lig_id)
            return Pose(100, 0, {}, 300)
        return self.poses[np.argmin([self.poses[i].rmsd for i in range(min(n, len(self.poses)))])]

class Docking:
    def __init__(self, dock_dir, ifp_dir, struct, chembl_ligs, load_fp, load_chembl):
        self.dock_dir = dock_dir
        self.ifp_dir = ifp_dir

        self.struct = struct
        self.ligands = {}
        self.load_docking(chembl_ligs, load_fp, load_chembl)

        self.mcss = {}

    def load_mcss(self, mcss_path, load_chembl):
        num_poses = {l:len(lig.poses) for l, lig in self.ligands.items()}
        self.mcss = parse_mcss(mcss_path, load_chembl, num_poses)

    def get_mcss_score(self, l1, p1, l2, p2):
        if (l1,l2) in self.mcss: return self.mcss[(l1,l2)][(p1,p2)]
        if (l2,l1) in self.mcss: return self.mcss[(l2,l1)][(p2,p1)]

        return None

    def load_docking(self, chembl_ligs, load_fp, load_chembl):
        ligands = []

        for pair in os.listdir(self.dock_dir):
            if pair[0] == '.': continue
            if pair[:6] == 'CHEMBL' and not load_chembl: continue
            l, s = pair.split('-to-')
            if 'self' == self.struct and s != l.split('_')[0]: continue
            elif s != self.struct: continue
            if not os.path.exists('{}/{}/{}_pv.maegz'.format(self.dock_dir, pair, pair)): continue
            ligands.append(l)

        for l in ligands:
            self.ligands[l] = Ligand(l)

            if l[:6] == 'CHEMBL':
                try:
                    self.ligands[l].ki = chembl_ligs[l].ki
                    self.ligands[l].stereo = chembl_ligs[l].valid_stereo
                except:
                    print 'no ki found for', l, 'drug?'
            else:
                self.ligands[l].ki = None
                self.ligands[l].stereo = None 

            st = self.struct
            if self.struct == 'self':
                st = l.split('_')[0]
            
            self.ligands[l].load_poses(self.dock_dir, self.ifp_dir, st, load_fp)

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
    def __init__(self, uniprot, data_dir):
        self.uniprot = uniprot

        self.pdb_ids = [l.split('_')[0] for l in os.listdir('{}/{}/ligands/unique'.format(data_dir, uniprot))]
        self.pdb_ids = sorted([l for l in self.pdb_ids
            if l+'_lig.mae' in os.listdir('{}/{}/structures/ligands'.format(data_dir, uniprot))])

        self.chembl_ligs = {}
        try:
            self.chembl_ligs = {'{}_lig'.format(c_id):c for c_id, c in load_chembl_proc('{}/{}'.format(data_dir, uniprot)).items()}
        except:
            print 'chembl not loaded'

        self.docking = {}
        self.true = {}

    def load_docking(self, dock_dir, ifp_dir, load_fp, load_chembl, structs):
        for st in structs:
            ld = Docking(dock_dir, ifp_dir, st, self.chembl_ligs, load_fp, load_chembl)
            if len(ld.ligands) > 0:
                self.docking[st] = ld

    def load_true(self, ifp_dir):
        ligs = [l.split('.')[0] for l in os.listdir(ifp_dir) if l[-2:] == 'fp']
        ligs = [l for l in ligs if len(l.split('-to-')) == 1 and len(l) > 0] 
        self.true = {l:Ligand(l) for l in ligs}
        for l, lig in self.true.items():
             lig.load_poses(None, ifp_dir, l)

class Dataset:
    def __init__(self, data_dir, prots):
        self.proteins = {}
        for u in prots:
            self.proteins[u] = Protein(u, data_dir)
        self.data_dir = data_dir

    def load_docking(self, gdir, fpdir, mdir, structs={},
                     load_fp=True, load_crystal=True, load_chembl=True, load_mcss=True):

        for u,p in self.proteins.items():
            dock_dir = '{}/{}/{}'.format(self.data_dir, u, gdir)
            ifp_dir = '{}/{}/{}'.format(self.data_dir, u, fpdir)
            st_list = p.pdb_ids
            if u in structs:
                st_list = [structs[u]]
            p.load_docking(dock_dir, ifp_dir, load_fp, load_chembl, st_list)

            if load_mcss:
                lig_feats = load_features('{}/{}'.format(self.data_dir, u))
                for s, st in p.docking.items():
                    mcss_dir = '{}/{}/{}/{}'.format(self.data_dir, u, mdir, s)
                    st.load_mcss(mcss_dir, load_chembl)

                    for l, lig in st.ligands.items():
                        lig.chrg = 'chrg' in lig_feats[l]
                        lig.ring = 'ring' in lig_feats[l]

            if load_crystal:
                p.load_true(ifp_dir, load_fp) 

    def assign_weights(self, w):
        for p_name, p in self.proteins.items():
            for d_name, d in p.docking.items():
                for l_name, l in d.ligands.items():
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

