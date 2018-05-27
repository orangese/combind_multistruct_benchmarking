import os
import sys

import numpy as np

sys.path.append('../3_analyze')
from containers import Dataset
from utils import export, show_side_by_side, show_features
from statistics import Statistics
from prob_opt import PredictStructs

class ScoreContainer:
    def __init__(self, data_dir, out_dir, prot, struct):
        self.prot = prot
        self.struct = struct
        self.root = '{}/{}/{}'.format(data_dir, prot, out_dir)
      
        sys.path.append(self.root)
        globals()['settings'] = __import__('settings')
 
        self.stats = self.init_stats()        

        self.predict_data = Dataset([prot], settings.data_dir, 
                                    settings.glide_dir, settings.ifp_dir, settings.mcss_dir)
        self.predict_data.proteins[prot].load([], struct, False, False, False)
        
        self.ps = PredictStructs(self.predict_data.proteins[prot].docking[struct], 
                                 self.stats.evidence, settings.features, 
                                 settings.num_poses, settings.t)
    
        self.results = {}

    def init_stats(self):
        stats_data = Dataset(settings.stats_prots, settings.data_dir, 
                             settings.glide_dir, settings.ifp_dir, settings.mcss_dir)
        
        stats_data.load({p:prot.lm.get_docked(pdb_only=True) 
                         for p,prot in stats_data.proteins.items()})

        return Statistics(stats_data, settings.stats_prots, 
                          settings.num_stats_ligs, settings.num_stats_poses, 
                          settings.features, settings.smooth, settings.normalize)

    def compute_results(self, q):
        prot = self.predict_data.proteins[self.prot]
        chembl_ligs = prot.lm.get_similar(q, settings.chembl_file, 
                                          num=settings.num_pred_chembl, 
                                          mcss_sort=settings.mcss_sort, struct=self.struct)
        self.predict_data.load({self.prot:[q]+chembl_ligs},{self.prot:[self.struct]}) 
        best_cluster, all_scores, all_rmsds = self.ps.max_posterior([q]+chembl_ligs, restart=15, sampling=3)
        result = self.ps.log_posterior(best_cluster)
 
        us_top = best_cluster[q]
    
        fpath = '{}/{}-to-{}.sc'.format(self.root, q, self.struct)
        with open(fpath,'a') as f:
            f.write('{},{}\n'.format(q, us_top))
            for lig, pose in sorted(best_cluster.items()):
                if lig == q: continue
                f.write('{},{}\n'.format(lig, pose))

            f.write('max_score,{}\n'.format(result))

    def load_results(self):
        self.results = {}
        for f in sorted(os.listdir(self.root)):
            if f.split('.')[-1] != 'sc' or f[0] == '.': continue
            q,s = f.split('.')[0].split('-to-')
            if s != self.struct: continue
            self.results[q] = {}
            with open('{}/{}'.format(self.root, f)) as outf:
                print f
                for line in outf:
                    line = line.strip().split(',')
                    if len(line) != 2: continue
                    if line[0] == 'max_score':
                        self.predict_data.load({self.prot:self.results[q].keys()},{self.prot:[self.struct]})
                        best_cluster = {l:p for l,p in self.results[q].items()}
                        print q, self.results[q].keys()
                        check_sc = self.ps.log_posterior(best_cluster)
                        if abs(check_sc - float(line[1])) > 0.0001:
                            print q, check_sc, float(line[1])
                        continue
                    try:
                        self.results[q][line[0]] = int(line[1])
                    except:
                        #errors.add(q)
                        print f
                        print line
    
    def get_rmsd(self, q, l, glide=False):
        if glide:
            return self.predict_data.proteins[self.prot].docking[self.struct].ligands[l].poses[0].rmsd
        return self.predict_data.proteins[self.prot].docking[self.struct].ligands[l].poses[self.results[q][l]].rmsd
    
    def show_fps(self, p, q, num_i=10, size=2):
        self.predict_data.assign_weights({2:1,3:1,4:1,6:1,11:0.005})
        l_list = [q]+sorted([l for l in self.results[q] if l!=q])
        us_top = self.ps.ligset.get_poses(self.results[q])
        glide_top = self.ps.ligset.get_poses({l:0 for l in l_list})
        show_side_by_side(us_top, glide_top, l_list, 
                          t1='Our Top Poses', t2='Glide Top Poses', num_i=num_i, size=size)
        
    def show_feature(self, p, q, f_name, show_prob=True, show_x=False, size=2):
        k,kdef = f_name, self.stats.features[f_name]
        l_list = [q]+sorted([l for l in self.results[p][q] if l!=q])
        
        x1, log_p1 = self.ps.x(self.results[p][q],k,kdef,lig_order=l_list)
        x2, log_p2 = self.ps.x({l:0 for l in l_list},k,kdef,lig_order=l_list)
        
        us_top = self.ps.ligset.get_poses(self.results[p][q])
        glide_top = self.ps.ligset.get_poses({l:0 for l in l_list})
        
        if show_prob and np.sum(log_p1) != 0:
            print k, 'probability matrix'
            minval = min(np.min(log_p1[np.nonzero(log_p1)]),np.min(log_p2[np.nonzero(log_p2)]))
            maxval = max(np.max(log_p1[np.nonzero(log_p1)]),np.max(log_p2[np.nonzero(log_p2)]))

            show_features(us_top, log_p1, glide_top, log_p2, l_list, 
                          'us top','glide top',size=size, mi=minval,ma=maxval)

        if show_x:
            print k, 'x_k matrix'
            minval = min(np.min(x1),np.min(x2))
            maxval = max(np.max(x1),np.max(x2))

            show_features(us_top, x1, glide_top, x2, l_list, 
                          'us top','glide top',size=size,mi=minval,ma=maxval)
            
    def export_poses(self, p, q):
        # this will show up in /scratch/PI/rondror/jbelk/method/outputs
        us_top = self.ps[p].ligset.get_poses(self.results[p][q])
        glide_top = self.ps[p].ligset.get_poses({l:0 for l in self.results[p][q]})
        export(self.data_dir, us_top, '{}_us_apr'.format(q), p, 
               struct=settings.struct_dict[p], verbose=False, glide_dir=settings.glide_dir)
        export(self.data_dir, glide_top, '{}_glide_apr'.format(q), p, 
               struct=settingsstruct_dict[p], verbose=False, glide_dir=settings.glide_dir)


