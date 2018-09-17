import numpy as np
from pairs import LigPair

class DensityEstimate:
    '''
    Computes and stores density estimates f(x) for input values x.
    
    DensityEstimates are stored in the form of a probability density
    to get a value in the form of counts at a position, just multiply
    by n_samples.

    DensityEstimates can be averaged and two can be combined to get a ratio.
    '''
    def __init__(self, points=100, domain = None, sd = 1.0, reflect = True,
                 n_samples = 0, out_of_bounds = None):
        '''
        points (int): number of values at which to compute density.
        domain ((float, float)): range of values at which to compute density
            if left as None, use min and max of input data.
        sd (float): standard deviation of gaussian kernel.
        reflect (bool): If True compute density for domain + left and right
            flanks, then reflect flanks and add to center. (This provides
            better behaviour at boundaries.)
        n_samples (int): (Weighted) number of samples used to compute density.
        out_of_bounds (None or float): If asked to return a density for a
            value outside the domain, if None, return closest density, else
            return this value.
        '''
        self.points = points
        self.out_of_bounds = out_of_bounds
        self.sd = sd
        self.n_samples = n_samples
        self.reflect = reflect
        if domain is not None:
            self.x = np.linspace(domain[0], domain[1], points)
            self.fx = np.zeros(self.x.shape)
        else:
            self.x = None
            self.fx = None

    # I/O
    # File format: first line -> n_samples, sd, reflect, remaining -> x, fx
    def __str__(self):
        return '\n'.join([','.join([str(self.n_samples), str(self.sd),
                                    str(self.reflect)])]
                         + ['{},{}'.format(_x, _fx)
                            for _x, _fx in zip(self.x, self.fx)])

    def write(self, fname):
        with open(fname, 'w') as fp:
            fp.write(str(self))

    @classmethod
    def read(cls, fname):
        x, fx = [], []
        with open(fname) as fp:
            n_samples, sd, reflect = fp.readline().strip().split(',')
            for line in fp:
                _x, _fx = line.strip().split(',')
                x += [float(_x)]
                fx += [float(_fx)]
        de = DensityEstimate(points = len(x),
                             n_samples=float(n_samples),
                             sd = float(sd),
                             reflect = bool(reflect))

        de.x, de.fx = np.array(x), np.array(fx)
        return de

    # Core methods
    def __call__(self, x):
        '''
        Returns f(x) for the given value of x by linear interpolation.
        If x is out of functions domain, return the closest response
        and print a warning if self.out_of_bounds is None or else
        self.out_of_bounds.
        '''
        if self.x is None: return 0
        if x < self.x[0]:
            print('Warning: value out of domain of function')
            if self.out_of_bounds is None: return self.fx[0]
            return self.out_of_bounds
        if x > self.x[-1]:
            print('Warning: value out of domain of function')
            if self.out_of_bounds is None: return self.fx[-1]
            return self.out_of_bounds

        idx = 0 # One to right of x
        while x > self.x[idx]: idx += 1 
        
        if not idx: return self.fx[0] # Avoid IndexError
        
        d = self.x[idx] - self.x[idx-1]
        return (   self.fx[idx-1]  * (self.x[idx] - x)/d
                 + self.fx[idx]    * (x - self.x[idx-1])/d)

    def _gauss(self, mean, x):
        '''
        Return PDF of N(mean, sd**2) at x.
        '''
        return (np.exp(-0.5*((x - mean)/self.sd)**2)
                / (self.sd*np.sqrt(2*np.pi)))

    def _kde(self, X, weights):
        '''
        Returns density estimate at each point in self.x for the input data X
        weighted by weights.
        '''
        self.fx = []
        for mean in self.x:
            kernel = np.array([self._gauss(mean, _X) for _X in X])
            self.fx += [(weights*kernel).sum()]
        self.fx = np.array(self.fx)


    def fit(self, X, weights = 1):
        '''
        Given an array of values X and weights weights,
        compute a density estimate with standard deviation self.sd.
        If reflect, compute densities for each flank and add
        computed densities back to the center.
        If hist, normalize so that area under the curve is equal to 
        '''
        if self.x is None:
           self.x = np.linspace(X.min(), X.max(), self.points)

        if self.reflect:
            r = self.x[-1]-self.x[0]
            self.x = np.hstack([self.x-r, self.x, self.x+r])

        self._kde(X, weights)
        
        if self.reflect:
            # left, center, right
            print self.x.shape, self.fx.shape
            self.fx = (  self.fx[self.points:0:-1]
                       + self.fx[self.points:2*self.points]
                       + self.fx[-1:2*self.points-1:-1])
            self.x = self.x[self.points:2*self.points]

        self.fx *= (self.x.shape[0] / (self.x[-1]-self.x[0])) / self.fx.sum()
        self.n_samples = (weights*np.ones(X.shape)).sum()
        return self

    # Combining instances
    def ratio(self, other, prob = True):
        '''
        Returns a new function representing the ratio of self
        over the other function. The domain of the new function
        covers the domain of both input functions with the same
        number of points as self. If prob is False weight ratio
        by number of samples.

        It would be interesting to evaluate if adding some shrinkage
        prior would give better results!
        '''
        de = DensityEstimate(points = self.points, out_of_bounds=1,
                             n_samples = self.n_samples)
        de.x = np.linspace(min(self.x[0], other.x[0]),
                           max(self.x[-1], other.x[-1]),
                           self.points)
        de.fx = np.array([self(x) / other(x) for x in de.x])
        if not prob: de.fx *= self.n_samples / float(other.n_samples)
        return de

    def average(self, other, weight = True):
        '''
        Returns a new function representing the average of
        the self and other functions. The domain of the new
        function covers the domain of both input functions with
        the same number of points as self. If weight, weight
        average by number of samples.
        '''
        assert self.reflect == other.reflect, "Either reflect or don't."
        if not other.n_samples:
            return self
        if not self.n_samples:
            return other
        
        de = DensityEstimate(points = self.points,
                             n_samples = self.n_samples + other.n_samples,
                             out_of_bounds=self.out_of_bounds,
                             reflect = self.reflect)
        de.x = np.linspace(min(self.x[0], other.x[0]),
                           max(self.x[-1], other.x[-1]),
                           self.points)

        if weight:
            de.fx = np.array([(  self(x) *  self.n_samples
                               + other(x)* other.n_samples)
                              / (self.n_samples+other.n_samples)
                              for x in de.x])
        else:
            de.fx = np.array([(self(x) + other(x)) / 2.0
                              for x in de.x])
        return de

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
