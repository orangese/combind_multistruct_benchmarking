import numpy as np
from pairs import LigPair

def average_distributions(d_list):
    all_x = set()
    all_raw = []
    for d in d_list:
        all_x.update(d.prob.keys())
        if d.raw is not None:
            all_raw += d.raw
        else:
            all_raw = None
    all_x = sorted(all_x)
    p_x = [np.mean([d.prob.get(x,0) for d in d_list]) for x in all_x]
    return Distribution(all_x, p_x, raw=all_raw)

def combine(dist, subdist):
    for i in dist:
        for k in dist[i]:
            all_d = [sv.dist[i][k] for sk,sv in subdist.items() if sv.dist[i][k] is not None]
            if len(all_d) == 0: continue
            dist[i][k] = average_distributions(all_d)

class Distribution:
    def __init__(self, x, p_x, raw=None):

        self.resolution = 2#int(np.log10(1/(x[1] - x[0])))
        self.prob = {round(x[i], self.resolution):p_x[i] for i in range(len(x))}
        self.dx = x[1] - x[0]
        self.raw = raw
        if type(raw) is dict:
            self.raw = []
            for ke, freq in raw.items(): 
                self.raw += freq*[ke]
 
    def evaluate(self, x):
        return self.prob[round(x,self.resolution)]

    def mean(self):
        return sum([x*px*self.dx for x,px in self.prob.items()])

class Statistics:
    def __init__(self, ligs, structs, k_list, normalize=True):
        self.proteins = {p:ProteinStatistics(p,structs[p],l,k_list) for p,l in ligs.items()}
        self.normalize = normalize
        self.k_list = k_list
        self.ind = [-1,0,1]
        self.dist = {i: {k: None for k in k_list} for i in self.ind}

    def create(self, dataset, samples, smooth, num_poses, hack=False):
        for p, pstats in self.proteins.items():
            prot = dataset.proteins[p]
            pstats.create(prot.docking[pstats.st], prot.lm.mcss, samples, smooth, num_poses, self.normalize, hack)
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

    def create(self, docking, mcss, samples, smooth, num_poses, normalize, hack):
        for (l1,l2), pstats in self.pairs.items():
            lp = LigPair(docking.ligands[l1], docking.ligands[l2],
                         self.k_list, mcss, num_poses, normalize)
            lp.init_pose_pairs()
            pstats.create(lp, samples, smooth, normalize, hack)
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
        self.x = {k:{i:{} for i in self.ind} for k in k_list}

    def create(self, lig_pair, samples=10**4, smooth=0.02, normalize=True, hack=False):
        for k in self.k_list:
            #x_k = {i: {} for i in self.ind}
            for (r1,r2),pp in lig_pair.pose_pairs.items():
                pp_x = lig_pair.get_feature(k, r1, r2)
                if pp_x is not None:
                    self.x[k][pp.correct()][pp_x] = self.x[k][pp.correct()].get(pp_x, 0) + 1#.append(pp_x)
                    self.x[k][-1][pp_x] = self.x[k][-1].get(pp_x, 0) + 1
                    if hack:
                        f = lambda x: sum([freq for ind,freq in x.items()])
                        if min(f(self.x[k][1]),f(self.x[k][0])) >= 4 and f(self.x[k][-1]) >= 250: break

            if len(self.x[k][-1]) == 0: continue
            if normalize: domain = np.linspace(-1,2,3*samples+1)    
            elif k == 'contact': domain = np.arange(-10,max(self.x[k][-1])+10,1)
            else: domain = np.arange(-0.5,max(self.x[k][-1])+0.5,0.01)

            var = smooth
            if not normalize:
                tot1 = float(sum([v for ke,v in self.x[k][-1].items()]))#/float(len(x_k[i]))
                mean = sum([ke*v/tot1 for ke,v in self.x[k][-1].items()])
                var = smooth*sum([((ke-mean)**2) * v/tot1 for ke,v in self.x[k][-1].items()])
                if var == 0:
                    print('variance error')
                    var = 0.1

            for i in self.ind:
                if len(self.x[k][i]) == 0: continue 
                tot = float(sum([v for ke,v in self.x[k][i].items()]))#/float(len(x_k[i]))
                gauss = lambda x,u: (1/np.sqrt(2*np.pi*var))*np.exp(-(x-u)**2 / (2*var))
                p_domain = [sum([gauss(j,v)*freq/tot for v,freq in self.x[k][i].items()]) for j in domain]

                self.dist[i][k] = Distribution(domain, p_domain, raw=self.x[k][i])
                #self.dist[i][k].raw = self.x[k][i]

    def read(self, data_dir, stats_dir):
        for k in self.k_list:
            out_f = '{}/{}/stats/{}/{}-{}-to-{}-{}.txt'.format(data_dir, 
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

    from shared_paths import shared_paths
    from containers import Dataset
    data = Dataset(shared_paths,[prot])
    data.load({prot:[l1,l2]}, load_fp=True, load_mcss=True)
    lm = data.proteins[prot].lm

    k_list = ['sb1','sb2','sb3','mcss','hbond','pipi','contact']
    alls = Statistics({prot:[l1,l2]}, {prot:lm.st}, k_list, normalize=True)
    alls.create(data,100,0.005,100)

    for k in k_list:
        out_f = '{}/{}/stats/{}/{}-{}-to-{}-{}.txt'.format(shared_paths['data'],prot,shared_paths['stats'],l1,l2,lm.st,k)    
        alls.write(out_f, k)
