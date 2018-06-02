import numpy as np

from pairs import LigPair

def average_distributions(d_list):
    all_x = set([])
    for d in d_list:
        all_x.update(d.prob.keys())
    all_x = sorted(list(all_x))
    p_x = [np.mean([d.prob.get(x,0) for d in d_list]) for x in all_x]
    return Distribution(all_x, p_x)

def combine(dist, subdist):
    for i in dist:
        for k in dist[i]:
            all_d = [sv.dist[i][k] for sk,sv in subdist.items() if sv.dist[i][k] is not None]
            if len(all_d) == 0: continue
            dist[i][k] = average_distributions(all_d)

class Distribution:
    def __init__(self, x, p_x):
        self.prob = {x[i]:p_x[i] for i in range(len(x))}
        self.resolution = int(np.log10(1/(x[1] - x[0])))

    def evaluate(self, x):
        return self.prob[round(x,self.resolution)]

class AllStatistics:
    def __init__(self, prots, features):
        self.proteins = {p:ProteinStatistics(l,features) for p,l in prots.items()}
        self.features = features
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in features} for i in self.ind}

    def create(self, dataset, samples, smooth, num_poses):
        for p, pstats in self.proteins.items():
            d, docking = dataset.proteins[p].docking.items()[0]
            pstats.create(docking, samples, smooth, num_poses)
        combine(self.dist, self.proteins)

    def write(self, data_dir, out_dir):
        for p, pstats in self.proteins.items():
            pstats.write('{}/{}/{}'.format(data_dir, p, out_dir))

    def read(self, data_dir, out_dir):
        for p, pstats in self.proteins.items():
            pstats.read('{}/{}/{}'.format(data_dir, p, out_dir))
        combine(self.dist, self.proteins)

class ProteinStatistics:
    def __init__(self, ligs, features):
        self.ligs = sorted(ligs)
        self.pairs = {}
        for i1,l1 in enumerate(self.ligs):
            for l2 in self.ligs[i1+1:]:
                self.pairs[(l1,l2)] = LigPairStatistics(l1,l2,features)

        self.features = features
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in features} for i in self.ind}

    def create(self, docking, samples, smooth, num_poses):
        for (l1,l2), pstats in self.pairs.items():
            lp = LigPair(docking.ligands[l1], docking.ligands[l2],
                         self.features, docking.mcss, num_poses, True)
            lp.init_pose_pairs()
            pstats.create(lp, samples, smooth)
        combine(self.dist, self.pairs)

    def write(self, out_dir):
        for (l1,l2), pstats in self.pairs.items():
            pstats.write(out_dir)

    def read(self, out_dir):
        for (l1,l2), pstats in self.pairs.items():
            pstats.read(out_dir)
        combine(self.dist, self.pairs)

class LigPairStatistics:
    def __init__(self, l1, l2, features):
        self.l1 = l1
        self.l2 = l2
        self.features = features
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in features} for i in self.ind}

    def create(self, lig_pair, samples=10**4, smooth=0.02):
        domain = np.linspace(-1,2,3*samples+1)
        for k in self.features:
            x_k = {i: [] for i in self.ind}
            for (r1,r2),pp in lig_pair.pose_pairs.items():
                pp_x = lig_pair.get_feature(k, r1, r2)
                if pp_x is not None:
                    x_k[pp.correct()].append(pp_x)
                    x_k[-1].append(pp_x)
            
            for i in self.ind:
                if len(x_k[i]) == 0: continue
                w = 1/float(len(x_k[i]))
                var = smooth
                gauss = lambda x,u: (w/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
                p_domain = [sum([gauss(j,v) for v in x_k[i]]) for j in domain]
                self.dist[i][k] = Distribution(domain, p_domain)

    def write(self, out_path):
        for k in self.features:
            out_f = '{}/{}-{}-{}.txt'.format(out_path, self.l1, self.l2, k)    
            with open(out_f, 'w') as f:
                for i in self.ind:
                    d = self.dist[i][k]
                    for x in sorted(d.prob.keys()):
                        f.write('{},{},{}\n'.format(i,x,d.prob[x]))

    def read(self, out_path):
        for k in self.features:
            out_f = '{}/{}-{}-{}.txt'.format(out_path, self.l1, self.l2, k)    
            temp_map = {i:{} for i in self.ind}
            with open(out_f, 'w') as f:
                for line in f:
                    i,x,p_x = line.strip().split(',')
                    temp_map[int(i)][float(x)] = float(p_x)
            for i in self.ind:
                x = sorted(temp_map[i].keys())
                p_x = [temp_map[i][xval] for xval in x]
                self.dist[i][k] = Distribution(x, p_x)


