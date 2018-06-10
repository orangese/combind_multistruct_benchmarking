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

class Statistics:
    def __init__(self, ligs, structs, k_list, normalize=True):
        self.proteins = {p:ProteinStatistics(p,structs[p],l,k_list) for p,l in ligs.items()}
        self.normalize = normalize
        self.k_list = k_list
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in k_list} for i in self.ind}

    def create(self, dataset, samples, smooth, num_poses):
        for p, pstats in self.proteins.items():
            prot = dataset.proteins[p]
            pstats.create(prot.docking[pstats.st], samples, smooth, num_poses, self.normalize)
        combine(self.dist, self.proteins)

    def read(self, data_dir, stats_dir):
        for p, pstats in self.proteins.items():
            self.proteins[p].read(data_dir, stats_dir)
        combine(self.dist, self.proteins)
 
    def write(self, out_f, k):
        with open(out_f, 'w') as f:
            for i in self.ind:
                d = self.dist[i][k]
                if d is None:
                    f.write('None\n')
                    break
                for x in sorted(d.prob.keys()):
                    f.write('{},{},{}\n'.format(i,x,d.prob[x]))

    def evaluate(self, k, x_k, i):
        if self.dist[i][k] is None:
            return 1
        return self.dist[i][k].evaluate(x_k)

class ProteinStatistics:
    def __init__(self, prot, st, ligs, k_list):
        self.prot = prot
        self.st = st
        self.ligs = sorted(ligs)
        self.pairs = {}
        for i1,l1 in enumerate(self.ligs):
            for l2 in self.ligs[i1+1:]:
                self.pairs[(l1,l2)] = LigPairStatistics(prot,st,l1,l2,k_list)

        self.k_list = k_list
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in k_list} for i in self.ind}

    def create(self, docking, samples, smooth, num_poses, normalize):
        for (l1,l2), pstats in self.pairs.items():
            lp = LigPair(docking.ligands[l1], docking.ligands[l2],
                         self.k_list, docking.mcss, num_poses, normalize)
            lp.init_pose_pairs()
            pstats.create(lp, samples, smooth, normalize)
        combine(self.dist, self.pairs)

    def read(self, data_dir, stats_dir):
        for (l1,l2), pstats in self.pairs.items():
            pstats.read(data_dir, stats_dir)
        combine(self.dist, self.pairs)

class LigPairStatistics:
    def __init__(self, prot, st, l1, l2, k_list):
        self.prot = prot
        self.st = st
        self.l1 = l1
        self.l2 = l2
        self.k_list = k_list
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in k_list} for i in self.ind}

    def create(self, lig_pair, samples=10**4, smooth=0.02, normalize=True):
        for k in self.k_list:
            x_k = {i: {} for i in self.ind}
            for (r1,r2),pp in lig_pair.pose_pairs.items():
                pp_x = lig_pair.get_feature(k, r1, r2)
                if pp_x is not None:
                    x_k[pp.correct()][pp_x] = x_k[pp.correct()].get(pp_x, 0) + 1#.append(pp_x)
                    x_k[-1][pp_x] = x_k[-1].get(pp_x, 0) + 1
            if normalize: domain = np.linspace(-1,2,3*samples+1)    
            else: domain = np.linspace(min(x_k[-1])-1,max(x_k[-1])+1,3*samples+1)    
            for i in self.ind:
                if len(x_k[i]) == 0: continue
                tot = float(sum([v for ke,v in x_k[i].items()]))#/float(len(x_k[i]))
                var = smooth
                gauss = lambda x,u: (1/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
                p_domain = [sum([gauss(j,v)*freq/tot for v,freq in x_k[i].items()]) for j in domain]
                self.dist[i][k] = Distribution(domain, p_domain)

    def read(self, data_dir, stats_dir):
        for k in self.k_list:
            out_f = '{}/{}/{}/{}-{}-to-{}-{}.txt'.format(data_dir, 
                self.prot, stats_dir, self.l1, self.l2, self.st, k)    
            temp = readf(out_f, k, self.ind)
            for i in self.ind:
                self.dist[i][k] = temp[i]

def readf(out_f, k, ind):
    temp_map = {i:{} for i in ind}
    with open(out_f) as f:
        for line in f:
            if line.strip() == 'None': 
                temp_map = {i:None for i in ind}
                break
            i,x,p_x = line.strip().split(',')
            temp_map[int(i)][float(x)] = float(p_x)
    tr = {}
    for i in ind:
        if temp_map[i] is None:
            tr[i] = None
            continue
        x = sorted(temp_map[i].keys())
        p_x = [temp_map[i][xval] for xval in x]
        tr[i] = Distribution(x, p_x)
    return tr

if __name__ == '__main__':
    import sys

    script_path, prot, l1, l2 = sys.argv
    code_path = '/'.join(script_path.split('/')[:-2])
    for i in ['1_dock','2_fp','3_analyze','4_score']:
        sys.path.append(code_path+'/'+i)

    from containers import Dataset, data_dir, stats_dir
    data = Dataset([prot])
    data.load({prot:[l1,l2]}, load_fp=True, load_mcss=True)
    lm = data.proteins[prot].lm

    k_list = ['sb1','sb2','sb3','mcss','hbond','pipi','contact']
    alls = Statistics({prot:[l1,l2]}, {prot:lm.st}, k_list)
    alls.create(data,10**2,0.02,100)

    for k in k_list:
        out_f = '{}/{}/{}/{}-{}-to-{}-{}.txt'.format(data_dir, prot, stats_dir, l1, l2, lm.st, k)    
        alls.write(out_f, k)



