import numpy as np
import math

from pairs import PosePair, LigPair
from stats_plots import stats_hist

class Statistics:
    def __init__(self, all_data, proteins, num_ligands, num_poses, features, smooth):
        self.all_data = all_data
        self.proteins = {p:{} for p in proteins} # maps proteins to ligpairs
        self.ligands = {p:prot.lm.pdb for p, prot in all_data.proteins.items()}
        self.features = features
        self.smooth = smooth
        self.evidence = self.create_distributions(num_ligands, num_poses)

    def create_distributions(self, num_ligs, num_poses):#_pairs):
        all_pairs = []
        for p in self.proteins:
            s,struct = self.all_data.proteins[p].docking.items()[0]
            #valid_ligs = self.ligands[p][:num_ligs]
            #for l in self.ligands[p]:
            #    valid_ligs[l] = [pose.rank for pose in struct.ligands[l].poses[:num_poses]]
            #    if len(valid_ligs) == num_ligs: break
            if len(self.ligands[p]) < num_ligs: 
                print 'warning, only {} ligands found for {}'.format(len(self.ligands[p]), p)
                num_ligs = len(self.ligands[p])

            for i1 in range(num_ligs):
                for i2 in range(i1+1, num_ligs):
                    l1,l2 = self.ligands[p][i1], self.ligands[p][i2]
                    lp = LigPair(struct.ligands[l1],struct.ligands[l2],
                                 self.features,struct.mcss, num_poses, num_poses)
                    self.proteins[p][(l1,l2)] = lp

        return Evidence(self.proteins, self.features, self.smooth)

    def show_stats(self, f_name, raw=True, smoothed=True, conditional=True):
        stats_hist(self.evidence, f_name, raw, smoothed, conditional)

    def show_stats_by_pair(self, f_name, raw=True, smoothed=True, conditional=True):
        for p,lps in self.proteins.items():
            for (l1,l2),lp in sorted(lps.items()):
                print l1, l2
                new_ev = Evidence({p:{(l1,l2):lp}}, {f_name:self.features[f_name]}, self.smooth)
                stats_hist(new_ev, f_name, raw, smoothed, conditional)

class Evidence:
    def __init__(self, all_data, features, smooth=0.25):
        self.all_data = all_data
        self.features = features
        self.smooth = smooth

        # 0: decoy
        # 1: correct
        # -1: all
        self.ind = [-1,0,1]
        self.raw_data = self.init_raw_data() 
        self.mean, self.std = self.init_mean_std()
        
        self.mult = {i:self.get_freq_map(i) for i in self.ind}
        self.cache = {i:{f:{} for f in features} for i in self.ind}

    def init_mean_std(self):
        mean = {f_name:{i:None for i in self.ind} for f_name in self.features}
        std = {f_name:{i:None for i in self.ind} for f_name in self.features}
        for f_name in self.features:
            if len(self.raw_data[f_name][-1]) == 0: continue
            for i in self.ind:
                mean[f_name][i] = np.mean(self.raw_data[f_name][i])
                std[f_name][i] = (np.mean([(val-mean[f_name][i])**2 for val in self.raw_data[f_name][i]]))**0.5
        return mean, std

    def get_freq_map(self, class_ind):
        freq_map = {f:{} for f in self.features}
        for f_name, f_def in self.features.items():
            for x_val in self.raw_data[f_name][class_ind]:
                freq_map[f_name][x_val] = freq_map[f_name].get(x_val, 0) + 1
            freq_map[f_name] = freq_map[f_name].items()
        return freq_map
        
    def init_raw_data(self):
        vals = {f_name:{i:[] for i in self.ind} for f_name in self.features}
        for f_name in self.features:
            for p, lps in self.all_data.items():
                for key,lp in lps.items():
                    for (r1,r2),pp in lp.pose_pairs.items():
                        pp_x = lp.get_feature(f_name, r1, r2)
                        if pp_x is not None:
                            vals[f_name][pp.correct()].append(pp_x)
                            vals[f_name][-1].append(pp_x)
        return vals
        
    def evaluate(self, f_name, f_val, pair_def):
        cache = self.cache[pair_def][f_name]
        freq_map = self.mult[pair_def][f_name]

        if f_val not in cache:
            w = 1/float(sum([mult for f,mult in freq_map]))
            var = self.smooth
            gauss = lambda x,u: (w/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
            cache[f_val] = sum([mult*gauss(f_val,f) for f,mult in freq_map])

        return cache[f_val]




