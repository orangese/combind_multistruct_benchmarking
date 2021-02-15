import numpy as np
from sklearn.linear_model import LinearRegression as _LinearRegression
from scipy.stats import norm
import matplotlib.pyplot as plt

class GaussianKernel:
    def __init__(self, sd=1.0, cutoff=float('inf')):
        self.sd = sd
        self.cutoff = cutoff

    def weights(self, X, center):
        w = (  (-self.cutoff*self.sd < X-center)
             & (X-center < self.cutoff*self.sd)).astype(float)
        w[w==1] = norm.pdf(X[w==1], loc=center, scale=self.sd)
        w /= 1 - 2*norm.cdf(-self.sd*self.cutoff, loc=0, scale=self.sd)
        return w

class LinearRegression(_LinearRegression):
    def zero_mean(self, X):
        y = self.predict(X)
        self.intercept_ -= y.mean()
        return self

    def inspect(self):
        print('Linear Regression: coeff={}, inter={}'.format(self.coef_[0], self.intercept_))

class FrozenLinearRegression(LinearRegression):
    def __init__(self, coeff):
        self.coef_ = np.array([coeff])
        self.intercept_ = 0

    def fit(self, X, y, w):
        return self

class LocalRegression:
    """
    Makes the given regressor local by applying the given kernel.
    Can be used either for regression or classification.
    """
    def __init__(self, kernel=GaussianKernel(), regressor=LinearRegression(),
                 points=100, domain=None):
        self.kernel = kernel
        self.regressor = regressor
        self.points = points
        self.domain = domain

        self.x = None
        self.fx = None
        self.n_samples = 0

    def fit(self, X, y, w=None):
        if w is None:
            w = np.ones(y.shape)

        assert len(X.shape) == 2, X.shape
        assert X.shape[1] == 1, X.shape
        assert len(y.shape) == 1, y.shape
        assert len(w.shape) == 1, w.shape
        assert X.shape[0] == y.shape[0] == w.shape[0]

        if self.domain is None:
            self.domain = (X.min(), X.max())

        self.x = np.linspace(self.domain[0], self.domain[1], self.points)
        self.fx = []
        for center in self.x:
            kernel = self.kernel.weights(X[:, 0], center)
            mask = kernel != 0
            r = self.regressor.fit(X[mask] - center,
                                   y[mask], (w*kernel)[mask])
            self.fx += [r.predict([[0]])[0]]
        self.fx = np.array(self.fx)
        self.n_samples = np.sum(w)
        return self

    def predict(self, X):
        assert len(X.shape) == 2, X.shape
        assert X.shape[1] == 1, X.shape
        return np.interp(X[:, 0], self.x, self.fx)

    def zero_mean(self, X):
        assert len(X.shape) == 2, X.shape
        assert X.shape[1] == 1, X.shape
        y = self.predict(X)
        self.fx -= y.mean()

    def inspect(self):
        print('Inspecting Local Regression:')
        plt.plot(self.x, self.fx)
        plt.show()

class LinearGAM:
    """
    Implements a linear generalized additive model

    f(X) = self.a + sum_i self.f_i(x_i)

    Model can have both linear and non-linear parts.
    """
    def __init__(self, smoother, convergence=0.00001, max_iterations=20):
        self.smoother = smoother
        self.convergence = convergence
        self.max_iterations = max_iterations

    def predict(self, X):
        y = np.zeros((X.shape[0],))
        y += self.a
        for col in X.columns:
            idx = ~np.isnan(X[col])
            y[idx] += self.f[col].predict(X[idx][[col]].to_numpy())
        return y

    def loss(self, X, y, w):
        return np.average((y - self.predict(X))**2, weights=w)

    def fit(self, X, y, w=None):
        if w is None:
            w = np.ones(y.shape)

        self.a = y.mean() # Intercept doesn't get updated.
        self.f = {col: LinearRegression().fit([[0], [1]], [0, 0])
                  for col in X.columns}

        loss = self.loss(X, y, w)
        for _ in range(self.max_iterations):
            print('Linear GAM MSE:', loss)
            for col in X.columns:
                # Compute residuals for current model with other columns
                other_cols = [c for c in X.columns if c != col]
                y_other = self.predict(X[other_cols])
                Y = y - y_other

                # Fit to residuals, ignoring places where there is no data.
                # Setting smoother to have zero mean gives stability.
                idx = ~np.isnan(X[col])
                self.f[col] = self.smoother(col).fit(X[idx][[col]].to_numpy(), Y[idx].to_numpy(), w[idx].to_numpy())
                #self.f[col].zero_mean((X[idx][[col]].to_numpy()))
            
            _loss = self.loss(X, y, w)
            if abs(loss - _loss) < self.convergence:
                break
            loss = _loss
        return self

    def inspect(self):
        print('Inspecting Linear GAM:')
        print('a={}'.format(self.a))
        for f, m in self.f.items():
            print(f)
            m.inspect()

class LogisticGAM:
    """
    Implements a logistic generalized additive model

    f(X) = 1 / (1 + exp(-self.a - sum_i self.f_i(x_i)))

    Model can have both linear and non-linear parts.
    """
    def __init__(self, smoother, convergence=0.00001, max_iterations=20,
                 regularization=0.01):
        self.smoother = smoother
        self.convergence = convergence
        self.max_iterations = max_iterations
        self.regularization = regularization

    def loss(self, X, y, w):
        p = self.predict_proba(X)
        return - np.average(y*np.log(p) + (1-y)*np.log(1-p), weights=w)

    def predict_proba(self, X):
        y = self.predict(X)
        return 1 / (1 + np.exp(-y))

    def predict(self, X):
        y = np.zeros((X.shape[0],))
        y += self.a
        for col in X.columns:
            idx = ~np.isnan(X[col])
            y[idx] += self.f[col].predict(X[idx][[col]].to_numpy())
        return y
        
    def fit(self, X, y, w=None):
        if w is None:
            w = np.ones(y.shape)

        self.a = np.log(y.mean() / (1 - y.mean()))
        self.f = {col: LinearRegression().fit([[0], [1]], [0, 0])
                  for col in X.columns}

        loss = self.loss(X, y, w)
        for _ in range(self.max_iterations):
            print('Logistic GAM cross-entropy loss:', loss)
            nu = self.predict(X)
            p = self.predict_proba(X)

            W = p * (1-p)
            W[W < self.regularization] = self.regularization

            z = nu + (y - p) / W

            lg = LinearGAM(self.smoother).fit(X, z, W*w)
            self.a = lg.a
            self.f = lg.f
            self.inspect()

            _loss = self.loss(X, y, w)
            if abs(loss - _loss) < self.convergence:
                break
            loss = _loss
        return self

    def inspect(self):
        print('Inspecting Logistic GAM:')
        print('a={}'.format(self.a))
        for f, m in self.f.items():
            print(f)
            m.inspect()
        print()
