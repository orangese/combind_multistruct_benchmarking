import numpy as np
import math

from pairs import PosePair, LigPair
from stats_plots import stats_hist

class Statistics:
    def __init__(self, all_data, proteins, num_ligands, num_pose_pairs, features):
        self.all_data = all_data
        self.proteins = {p:[] for p in proteins}
        self.features = features
        self.lig_pairs = {}
        self.evidence = self.create_distributions(num_ligands, num_pose_pairs)

    def create_distributions(self, num_ligs, num_pose_pairs):
        nnative = int(math.ceil(num_pose_pairs**0.5))
        all_pairs = []
        for p in self.proteins:
            prot = self.all_data.proteins[p]
            assert len(prot.docking) == 1 # should only be one structure
            for s, struct in prot.docking.items():
                valid_ligs = {}
                for l in prot.pdb_ids:
                    # must be cross docked
                    #if l == s: continue
                    l = '{}_lig'.format(l)
                    if l not in struct.ligands: continue
                    # and have some native/non native poses
                    native = [pose for pose in struct.ligands[l].poses[:100] if pose.rmsd <= 2]
                    decoy  = [pose for pose in struct.ligands[l].poses[:100] if pose.rmsd >  2]
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
                                     valid_ligs[l1],valid_ligs[l2],self.features,struct.mcss)
                        all_pairs.extend(lp.all_pose_pairs())
                        self.lig_pairs[(l1,l2)] = lp

        return Evidence(all_pairs, self.features)

    def show_stats(self, f_name, raw=True, smoothed=True):
        stats_hist(self.evidence, f_name, raw, smoothed)

    def show_stats_by_pair(self, f_name, raw=True, smoothed=True):
        for (l1,l2),lp in sorted(self.lig_pairs.items()):
            print l1, l2
            new_ev = Evidence(lp.all_pose_pairs(), {f_name:self.features[f_name]})
            stats_hist(new_ev, f_name, raw, smoothed)

class Evidence:
    def __init__(self, all_pairs, features):
        self.native = [p for p in all_pairs if p.native()]
        self.decoy = [p for p in all_pairs if not p.native()]
        
        self.features = features
        
        self.mean, self.std = self.init_mean_std()
        
        self.native_mult = self.get_freq_map('native')
        self.decoy_mult = self.get_freq_map('decoy')

        self.native_cache = {f:{} for f in features}
        self.decoy_cache = {f:{} for f in features}

    def init_mean_std(self):
        mean = {}
        std = {}
        for f_name,f_def in self.features.items():
            mean[f_name] = np.mean(self.raw_data(f_name))
            std[f_name] = (np.mean([(val-mean[f_name])**2 for val in self.raw_data(f_name)]))**0.5
        return mean, std

    def get_freq_map(self, pair_def):
        freq_map = {}
        for f_name, f_def in self.features.items():
            freq_map[f_name] = {}
            for x_val in self.raw_data(f_name, pair_def):
                freq_map[f_name][x_val] = freq_map[f_name].get(x_val, 0) + 1
            freq_map[f_name] = freq_map[f_name].items()
        return freq_map
        
    def raw_data(self, f_name, pair_def='native'):
        if pair_def == 'native':
            pairs = self.native
        elif pair_def == 'decoy':
            pairs = self.decoy
        else:
            pairs = self.decoy + self.native
            
        vals = [p.get_feature(f_name,self.features[f_name]) for p in pairs]
        return [v for v in vals if v is not None]
        
    def evaluate(self, f_name, f_val, pair_def='native', smooth=0.5):
        if pair_def == 'native':
            cache = self.native_cache[f_name]
            pairs = self.native
            freq_map = self.native_mult[f_name]
        elif pair_def == 'decoy':
            cache = self.decoy_cache[f_name]
            pairs = self.decoy
            freq_map = self.decoy_mult[f_name]

        if f_val not in cache:
            w = 1/float(sum([mult for f,mult in freq_map]))
            var = (smooth*self.std[f_name])**2
            gauss = lambda x,u: (w/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
            cache[f_val] = sum([mult*gauss(f_val,f) for f,mult in freq_map])

        return cache[f_val]




