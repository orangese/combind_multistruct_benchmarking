import os
import sys

import numpy as np

from containers import Dataset, LigandManager
from utils import export, show_side_by_side, show_features
from statistics import Statistics, readf
from prob_opt import PredictStructs

def read_settings(pth):
    tr = {}
    with open('{}/settings.py'.format(pth)) as f:
        for line in f:
            var,val = line.split('=')
            tr[var] = eval(val)
    return tr

class ScoreContainer:
    def __init__(self, root, prot, struct):
        self.prot = prot
        self.struct = struct
        self.root = root
      
        self.sett = read_settings(self.root)
        self.sp = self.sett['shared_paths']
 
        self.stats = self.init_stats()
        self.predict_data = Dataset(self.sp, [prot], {prot:struct})
        self.ps = PredictStructs(self.predict_data.proteins[prot], self.stats, 
            self.sett['k_list'], self.sett['num_poses'], self.sett['t'], self.sett['score_mode'])
    
        self.results = {}
        self.validate = {}
        self.details = set([])

    def init_stats(self):
        if os.path.exists('{}/stats_{}.txt'.format(self.root, self.sett['k_list'][0])):
            stats = Statistics({},{},self.sett['k_list'])
            for k in self.sett['k_list']:
                tmp = readf('{}/stats_{}.txt'.format(self.root, k), k, stats.ind)
                for i in stats.ind:
                    stats.dist[i][k] = tmp[i]
        else:
            stats_ligs = {}
            stats_st = {}
            for p in self.sett['stats_prots']:
                lm = LigandManager(self.sp,p)
                stats_ligs[p] = lm.docked(lm.pdb)[:self.sett['num_stats_ligs']]
                stats_st[p] = lm.st
            stats = Statistics(stats_ligs, stats_st, self.sett['k_list'])
            stats.read(self.sp['data'], self.sp['stats'])
        return stats

    def compute_results(self, q):
        if q not in self.results:
            prot = self.predict_data.proteins[self.prot]
            chembl_ligs = prot.lm.get_similar(q, self.sett['chembl_file'], 
                num=self.sett['num_pred_chembl'], struct=self.struct)
            self.predict_data.load({self.prot:[q]+chembl_ligs},{self.prot:[self.struct]}) 
            best_cluster, all_scores, all_rmsds = self.ps.max_posterior([q]+chembl_ligs, restart=15, sampling=3)
            self.results[q] = best_cluster
        return self.results[q]

    def load_results(self):
        self.results = {}
        for f in sorted(os.listdir(self.root)):
            if f.split('.')[-1] != 'sc' or f[0] == '.': continue
            q,s = f.split('.')[0].split('-to-')
            if s != self.struct: continue
            self.results[q] = {}
            with open('{}/{}'.format(self.root, f)) as outf:
                for line in outf:
                    line = line.strip().split(',')
                    if len(line) != 2: continue
                    if line[0] == 'max_score':
                        self.predict_data.load({self.prot:self.results[q].keys()},
                            {self.prot:[self.struct]}, load_fp=False, load_mcss=False)
                        self.validate[q] = float(line[1])
                        continue
                    try:
                        self.results[q][line[0]] = int(line[1])
                    except:
                        print f
                        print line

    def load_details(self, q):
        if q not in self.details:
            self.predict_data.load({self.prot:self.results[q].keys()},{self.prot:[self.struct]})
            self.details.add(q)    

    def validate_results(self,q):
        self.load_details(q)
        check_sc = self.ps.log_posterior(self.results[q])
        if abs(check_sc - self.validate[q]) > 0.0001:
            print check_sc, self.validate[q]
            return False
        return True

    def get_rmsd(self, q, l, glide=False):
        if glide:
            return self.predict_data.proteins[self.prot].docking[self.struct].ligands[l].poses[0].rmsd
        return self.predict_data.proteins[self.prot].docking[self.struct].ligands[l].poses[self.results[q][l]].rmsd
    
    def show_fps(self, q, num_i=10, size=2, w={2:1,3:1,4:1,6:1,11:0.005}):
        self.load_details(q)
        self.predict_data.assign_weights(w)
        l_list = [q]+sorted([l for l in self.results[q] if l!=q])
        us_top = self.ps.get_poses(self.results[q])
        glide_top = self.ps.get_poses({l:0 for l in l_list})
        show_side_by_side(us_top, glide_top, l_list, 
                          t1='Our Top Poses', t2='Glide Top Poses', num_i=num_i, size=size)
        
    def show_feature(self, q, k, show_prob=True, show_x=False, size=2, best=False):
        self.load_details(q)
        l_list = [q]+sorted([l for l in self.results[q].keys() if l != q])
        
        x1, log_p1 = self.ps.likelihood_and_feature_matrix(self.results[q],k,lig_order=l_list)
        x2, log_p2 = self.ps.likelihood_and_feature_matrix({l:0 for l in l_list},k,lig_order=l_list)
        
        title = 'glide top'
        us_top = self.ps.get_poses(self.results[q])
        glide_top = self.ps.get_poses({l:0 for l in l_list})

        if best:
            min_rmsd = np.argmin([p.rmsd for p in self.predict_data.proteins[self.prot].docking[self.struct].ligands[q].poses])
            glide_top = self.ps.get_poses(self.results[q])
            glide_top[q] = self.predict_data.proteins[self.prot].docking[self.struct].ligands[q].poses[p]
            title = 'min query rmsd'

        if show_prob and np.sum(log_p1) != 0:
            print k, 'probability matrix'
            minval = min(np.min(log_p1[np.nonzero(log_p1)]),np.min(log_p2[np.nonzero(log_p2)]))
            maxval = max(np.max(log_p1[np.nonzero(log_p1)]),np.max(log_p2[np.nonzero(log_p2)]))

            show_features(us_top, log_p1, glide_top, log_p2, l_list, 
                          'us top',title,size=size, mi=minval,ma=maxval)

        if show_x:
            print k, 'x_k matrix'
            minval = min(np.min(x1),np.min(x2))
            maxval = max(np.max(x1),np.max(x2))

            show_features(us_top, x1, glide_top, x2, l_list, 
                          'us top',title,size=size,mi=minval,ma=maxval)
            
    def export_poses(self, q):
        # this will show up in /scratch/PI/rondror/jbelk/method/outputs
        us_top = self.ps.get_poses(self.results[q])
        glide_top = self.ps.get_poses({l:0 for l in self.results[q]})
        export(self.sp['data'], us_top, '{}_us'.format(q), self.prot, 
               struct=self.struct, verbose=False, glide_dir=self.sp['docking'])
        export(self.sp['data'], glide_top, '{}_glide'.format(q), self.prot, 
               struct=self.struct, verbose=False, glide_dir=self.sp['docking'])

if __name__ == '__main__':
    script_path, s, p = sys.argv[:3]
    q_list = sys.argv[3:]

    sc = ScoreContainer(os.getcwd(), p, s)

    for k in sc.stats.k_list:
        out_f = '{}/stats_{}.txt'.format(sc.root, k)
        sc.stats.write(out_f, k)

    for q in q_list:
        best_cluster = sc.compute_results(q)
    
        fpath = '{}/{}-to-{}.sc'.format(sc.root, q, s)
        with open(fpath,'w') as f:
            f.write('{},{}\n'.format(q, best_cluster[q]))
            for lig, pose in sorted(best_cluster.items()):
                if lig == q: continue
                f.write('{},{}\n'.format(lig, pose))

            f.write('max_score,{}\n'.format(sc.ps.log_posterior(best_cluster)))


