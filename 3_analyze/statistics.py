import numpy as np
import math

from pairs import PosePair, LigPair
from stats_plots import stats_hist

class Statistics:
    def __init__(self, all_data, proteins, num_ligands, num_pose_pairs, features, smooth):
        self.all_data = all_data
        self.proteins = {p:[] for p in proteins}
        self.ligands = {p:prot.lm.get_pdb() for p, prot in all_data.proteins.items()}
        self.features = features
        self.lig_pairs = {}
        self.smooth = smooth
        self.evidence = self.create_distributions(num_ligands, num_pose_pairs)

    def create_distributions(self, num_ligs, num_pose_pairs):
        nnative = int(math.ceil(num_pose_pairs**0.5))
        all_pairs = []
        for p in self.proteins:
            struct = self.all_data.proteins[p].docking
            valid_ligs = {}
            for l in self.ligands[p]:
                native = [pose.rank for pose in struct.ligands[l].poses[:100] if pose.rmsd <= 2]
                decoy  = [pose.rank for pose in struct.ligands[l].poses[:100] if pose.rmsd >  2]
                if len(native) < nnative or len(decoy) < nnative: continue
                valid_ligs[l] = native[:nnative] + decoy[:nnative]
                self.proteins[p].append(l)
                if len(valid_ligs) == num_ligs: break
            else: 
                print 'warning, only {} ligands found for {}'.format(len(valid_ligs), p)
                num_ligs = len(valid_ligs)

            for i1 in range(num_ligs):
                for i2 in range(i1+1, num_ligs):
                    l1,l2 = self.proteins[p][i1], self.proteins[p][i2]
                    lp = LigPair(struct.ligands[l1],struct.ligands[l2],
                                 self.features,struct.mcss)
                    all_pairs.extend(lp.all_pose_pairs(valid_ligs[l1],valid_ligs[l2]))
                    self.lig_pairs[(l1,l2)] = lp

        return Evidence(all_pairs, self.features, self.smooth)

    def show_stats(self, f_name, raw=True, smoothed=True):
        stats_hist(self.evidence, f_name, raw, smoothed)

    def show_stats_by_pair(self, f_name, raw=True, smoothed=True):
        for (l1,l2),lp in sorted(self.lig_pairs.items()):
            print l1, l2
            new_ev = Evidence(lp.all_pose_pairs(), {f_name:self.features[f_name]})
            stats_hist(new_ev, f_name, raw, smoothed)

class Evidence:
    def __init__(self, all_pairs, features, smooth=0.25):
        self.all_pairs = all_pairs
        self.correct = {i:[p for p in all_pairs if p.correct() == i] for i in range(2)}

        self.smooth = smooth
        self.features = features
        
        self.mean, self.std = self.init_mean_std()
        
        self.mult = {i:self.get_freq_map(i) for i in range(2)}
        self.cache = {i:{f:{} for f in features} for i in range(2)}

    def init_mean_std(self):
        mean = {}
        std = {}
        for i in range(-1,2):
            mean[i] = {}
            std[i] = {}
            for f_name,f_def in self.features.items():
                mean[i][f_name] = np.mean(self.raw_data(f_name, i))
                std[i][f_name] = (np.mean([(val-mean[i][f_name])**2 for val in self.raw_data(f_name, i)]))**0.5
        return mean, std

    def get_freq_map(self, pair_def):
        freq_map = {}
        for f_name, f_def in self.features.items():
            freq_map[f_name] = {}
            for x_val in self.raw_data(f_name, pair_def):
                freq_map[f_name][x_val] = freq_map[f_name].get(x_val, 0) + 1
            freq_map[f_name] = freq_map[f_name].items()
        return freq_map
        
    def raw_data(self, f_name, pair_def=-1):
        pairs = self.correct.get(pair_def, self.all_pairs)   
        vals = [p.get_feature(f_name,self.features[f_name]) for p in pairs]
        return [v for v in vals if v is not None]
        
    def evaluate(self, f_name, f_val, pair_def):
        cache = self.cache[pair_def][f_name]
        pairs = self.correct[pair_def]
        freq_map = self.mult[pair_def][f_name]

        if f_val not in cache:
            w = 1/float(sum([mult for f,mult in freq_map]))
            var = (self.smooth*self.std[-1][f_name])**2
            gauss = lambda x,u: (w/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
            cache[f_val] = sum([mult*gauss(f_val,f) for f,mult in freq_map])

        return cache[f_val]




