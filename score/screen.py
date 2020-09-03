import numpy as np

def load_features_screen(features, gscore_fname, ifp_fname,
                         mcss_fname=None, shape_fname=None):
    single = np.load(gscore_fname)

    raw = {}
    for feature in features:
        if feature == 'mcss':
            raw['mcss'] = np.load(mcss_fname)
        elif feature == 'shape':
            raw['shape'] = np.load(shape_fname)
        else:
            raw[feature] = np.load(ifp_fname.format(feature))
    return single, raw

def screen(single, raw, stats, alpha):
    energies = {}
    for feature in raw:
        _raw = raw[feature]
        _stats = stats[feature]
        energies[feature] = (  np.log(_stats['native'](_raw))
                             - np.log(_stats['reference'](_raw)))

    pair_energy = 0
    for feature, energy in energies.items():
        pair_energy += energy

    pair_energy = pair_energy.mean(axis=1)
    combind_energy = pair_energy - alpha*single
    return combind_energy
