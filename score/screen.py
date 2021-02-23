import numpy as np
from utils import np_load

def load_features_screen(features, gscore_fname, ifp_fname,
                         mcss_fname=None, shape_fname=None):
    single = np.load(gscore_fname)

    raw = {}
    for feature in features:
        if feature == 'mcss':
            raw['mcss'] = np_load(mcss_fname)
        elif feature == 'shape':
            raw['shape'] = np_load(shape_fname)
        else:
            raw[feature] = np_load(ifp_fname.format(feature))
    return single, raw

def screen(single, raw, stats, alpha, weights=None):
    energies = {}
    for feature in raw:
        _raw = raw[feature]
        _stats = stats[feature]
        energies[feature] = (  np.log(_stats['native'](_raw))
                             - np.log(_stats['reference'](_raw)))

    pair_energy = 0
    for feature, energy in energies.items():
        pair_energy += energy

    n = pair_energy.shape[1]

    if weights is None:
        weights = np.ones(n)

    alpha /= 0.5 * n / (1 + (n-1)*0.5)

    pair_energy = (pair_energy*weights.reshape(1, -1)).mean(axis=1)
    combind_energy = pair_energy - alpha*single
    return combind_energy
