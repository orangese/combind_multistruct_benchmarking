import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
import sys

sys.path.append('/Users/jpaggi/Documents/combind/combind')
from shared_paths import proteins, shared_paths

pd.set_option("display.max_rows", 75)

families = {
    'GPCR': ['5HT2B', 'A2AR', 'B1AR', 'B2AR', 'CHRM3','SMO', 'MGLUR5'],
    'Kinase': ['BRAF', 'CDK2', 'CHK1', 'JAK2', 'PLK1', 'MAPK14', 'MEK1'],
    #'Ion Channel': ['TRPV1'],
    'Transporter': ['SLC6A4', 'GLUT1', 'DAT'],
    'Nuclear Receptor': ['NR3C2', 'NR3C1', 'AR', 'VDR', 'ERA'],
    'Peptidase': ['F2', 'F10', 'F11', 'PLAU', 'P00760', 'BACE1'],
    'Other': ['PYGM', 'PTPN1', 'BRD4', 'HSP90AA1', 'PDE10A', 'SIGMAR1', 'ELANE', 'TRPV1', 'DHFR']
}

drugs = {'GPCR': 0.33,
         'Kinase': 0.03,
         #'Ion Channel': 0.18,
         'Nuclear Receptor': 0.16,
         'Other': 0.20+0.18,
         'Peptidase': 0.03,
         'Transporter': 0.07}

family = pd.api.types.CategoricalDtype(['GPCR', 'Nuclear Receptor', 'Transporter', 'Peptidase', 'Kinase', 'Other'],
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

def add_correct(data, thresh = 2.0):
    data = data.copy()
    for col in data:
        if 'rmsd' in col:
            data[col.replace('rmsd', 'correct')] = data[col] < thresh
    return data

def filter_to_ubiquitous_ligands(data):
    level = [i for i, name in enumerate(data.index.names) if name not in ['family', 'protein', 'ligand']]
    ligands = set.intersection(*[set(group.index.get_level_values(-1).to_numpy())
                                 for name, group in data.groupby(level=level)])
    return data[data.index.get_level_values(-1).isin(ligands)]

def load(version, helpers, mcss):
    fnames = ('{}/{}/scores/{}/summary/{}.tsv'.format(shared_paths['data'], protein, version, helpers)
             for protein in proteins)
    data = pd.concat(pd.read_csv(fname, sep='\t')
                     for fname in fnames
                     if os.path.exists(fname))

    data['version'] = version
    data['helpers'] = helpers

    # Add mcss overlap with crystal ligand.
    with open(mcss, 'rb') as fp:
        mcss = pd.DataFrame(pickle.load(fp).items(), columns=['ligand', 'mcss']).set_index('ligand')
        data = data.join(mcss, on='ligand')
        
    # Add family
    reverse = {v:k for k, vs in families.items() for v in vs}
    data['family'] = data['protein'].apply(lambda x: reverse[x]).astype(family)

    data = data.set_index(['version', 'helpers', 'mode', 'features', 'alpha', 'n_ligs',
                           'family', 'protein', 'ligand']).sort_index()

    return data