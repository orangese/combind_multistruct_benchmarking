import numpy as np

class DensityEstimate:
    '''
    Computes and stores density estimates f(x) for input values x.
    
    DensityEstimates are stored in the form of a probability density
    to get a value in the form of counts at a position, just multiply
    by n_samples.

    DensityEstimates can be averaged and two can be combined to get a ratio.

    - Should explore ways to bias the ratio to 1 where there is little data.
    - Could consider collapsing values using dictionary when computing density
      to improve performance when there are many repeated items.
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
            kernel = self._gauss(mean, X)
            self.fx += [(weights*kernel).sum()]
        self.fx = np.array(self.fx)

    def uniform(self):
        self.n_samples = 0
        self.fx = np.ones(self.x.shape)
        self.fx /= self.x[-1]-self.x[0]
        return self

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

        if not X.shape[0]: return self.uniform()

        if self.reflect:
            if X.max() > self.x[-1] or X.min() < self.x[0]:
                print('Warning: Data out of domain of density estimate'
                      ' with reflected boundary conditions. Truncating'
                      ' data to be on specified domain.')
                X = X[(X <= self.x[-1])*(X>=self.x[0])]
                if not X.shape[0]: return self.uniform()
            r = self.x[-1] - self.x[0]
            self.x = np.hstack([self.x-r, self.x, self.x+r])

        self._kde(X, weights)
        
        if self.reflect:
            # left, center, right
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
        Returns a new function representing the ratio of self over the other
        function. The domain of the new function covers the domain of both
        input functions with the same number of points as self. If not prob 
        , weight ratio by number of samples.

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

    def _average(self, other):
        '''
        Returns a new function representing the average of the self and other
        functions. The domain of the new function covers the domain of both
        input functions with the same number of points as self.
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
                           max(self.x[-1], other.x[-1]), self.points)

        self_fx  = np.array([ self(x) for x in de.x])
        other_fx = np.array([other(x) for x in de.x])
        de.fx = ((self_fx * self.n_samples + other_fx * other.n_samples)
                 / (self.n_samples+other.n_samples))
        return de

    def copy(self):
        de = DensityEstimate()
        de.points = self.points
        de.out_of_bounds = self.out_of_bounds
        de.sd = self.sd
        de.n_samples = self.n_samples
        de.reflect = self.reflect
        de.domain = self.domain
        de.x = np.copy(self.x)
        de.fx = np.copy(self.fx)

    @classmethod
    def merge(cls, des, weight_equally = True):
        if weight_equally:
            for de in des:
                if de.n_samples:
                    de.n_samples = 1
        out = des[0]
        for de in des[1:]:
            out = out._average(de)
        return out
