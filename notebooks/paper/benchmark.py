import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
import sys
from glob import glob

pd.set_option("display.max_rows", 75)

families = {
    'GPCR': ['5HT2B', 'A2AR', 'B1AR', 'B2AR', 'SMO', 'MGLUR5'],
    'Ion Channel': ['P19491', 'P22756', 'Q05586-Q12879'],
    'Transporter': ['SLC6A4', 'GLUT1', 'DAT'],
    'Nuclear Receptor': ['NR3C2', 'NR3C1', 'AR', 'VDR', 'ERA'],
    'Peptidase': ['F2', 'F10', 'F11', 'PLAU', 'P00760', 'BACE1'],
    'Other': ['CDK2', 'PYGM', 'PTPN1', 'BRD4', 'HSP90AA1', 'PDE10A', 'SIGMAR1', 'ELANE', 'DHFR']
}

drugs = {'GPCR': 0.33,
         'Ion Channel': 0.18,
         'Nuclear Receptor': 0.16,
         'Other': 0.20+0.03,
         'Peptidase': 0.03,
         'Transporter': 0.07}

family = pd.api.types.CategoricalDtype(['GPCR', 'Ion Channel', 'Nuclear Receptor', 'Transporter',
                                        'Peptidase', 'Other'],
                                        ordered = True)

def dumbell_plot(data, metric, level = -2):
    label = data.index.names[level]
    f, ax = plt.subplots(figsize = (3, len(data)/3.1))
    data = data.sort_values([label, 'glide_{}'.format(metric)], ascending=[False, True])
    
    last = None
    i = -1
    yticks = []
    for label, values in data.iterrows():
        if level != -1 and last != label[level]:
            i += 2
        last = label[level]
        yticks += [i]
        
        g, c, b = values[['glide_{}'.format(metric),
                          'combind_{}'.format(metric),
                          'best_{}'.format(metric)]]
        plt.plot([g, c], [i, i], c = 'grey')
        plt.scatter(c, i, c = 'g')
        plt.scatter(g, i, c = 'm')
        i += 1
    
    yticklabels = list(data.index.get_level_values(-1))
    plt.yticks(yticks, yticklabels, fontsize = 12)
    if metric == 'correct':
        plt.xlim(0, 1)

def drug_average(family):
    assert 'protein' not in family.index.names
    assert 'ligand' not in family.index.names
    level = [i for i, name in enumerate(family.index.names) if name not in ['family']]
    weights = np.array([drugs[family] for family in family.index.get_level_values('family')])
    weighted = family * weights.reshape(-1, 1)
    return weighted.groupby(level=level).sum()

################################################################################

def add_correct(data, thresh = 2.0):
    data = data.copy()
    for col in data:
        if 'rmsd' in col:
            data[col.replace('rmsd', 'correct')] = round(data[col], 1) <= thresh
    return data

def filter_to_ubiquitous_ligands(data):
    level = [i for i, name in enumerate(data.index.names) if name not in ['family', 'protein', 'ligand']]
    ligands = set.intersection(*[set(group.index.get_level_values(-1).to_numpy())
                                 for name, group in data.groupby(level=level)])
    return data[data.index.get_level_values(-1).isin(ligands)]

def get_data(helpers, version, best=True):
    roots = ['bpp_data', 'ionchannels']
    data = pd.concat(load(version, helpers, root, '{}_mcss.pkl'.format(root)) for root in roots)
    data = add_correct(data, thresh = 2.0)
    data = data[data.mcss < 0.5]
    if best: data = data[data.best_correct]
    return data

def results(data, helpers, alpha=1.0, method='standard',
            features='mcss_contact_hbond_sb', ligs=20, aggregate='ligand'):
    ligand = data.xs((helpers, method, ligs, features, alpha),
                      level=('helpers', 'mode', 'n_ligs', 'features', 'alpha')).filter(regex='rmsd|correct')
    target = ligand.groupby(level=list(range(len(ligand.index.levels)-1))).mean().dropna()
    print(ligand.shape[0])
    if aggregate == 'target':
        family = target.groupby(level=list(range(len(target.index.levels)-1))).mean()
    elif aggregate == 'ligand':
        family = ligand.groupby(level=list(range(len(ligand.index.levels)-2))).mean()
    perf = drug_average(family)
    return perf['glide_correct'][0], perf['combind_correct'][0]

def load(version, helpers, root, mcss):
    data_root = '/Users/jpaggi/sherlock/oak/users/jpaggi'
    if root == 'bpp_data':
        version = 'stats41'
    if root == 'ionchannels':
        version = 'stats104'
    fnames = glob('{}/{}/*/scores/{}/summary/{}.tsv'.format(data_root, root, version, helpers))
    version = 'stats104'
    data = pd.concat(pd.read_csv(fname, sep='\t')
                     for fname in fnames)
    data = data[~data['protein'].isin(['O35433', 'P13569', 'O14649', 'TRPV1'])]

    data['version'] = version
    data['helpers'] = helpers

    # Add mcss overlap with crystal ligand.
    with open('{}/{}'.format(data_root, mcss), 'rb') as fp:
        mcss = pd.DataFrame(pickle.load(fp).items(), columns=['ligand', 'mcss']).set_index('ligand')
        data = data.join(mcss, on='ligand')
 
    # Add family
    reverse = {v:k for k, vs in families.items() for v in vs}
    data['family'] = data['protein'].apply(lambda x: reverse[x]).astype(family)

    data = data.set_index(['version', 'helpers', 'mode', 'features', 'alpha', 'n_ligs',
                           'family', 'protein', 'ligand']).sort_index()
    data.loc[data.index.get_level_values('protein') == 'TRPV1', 'mcss'] = 0
    return data