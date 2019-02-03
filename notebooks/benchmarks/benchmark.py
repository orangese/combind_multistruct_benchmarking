"""
Functions to plot ComBind v. Glide comparisons.
"""

import math
import numpy as np
from scipy.stats import ttest_rel
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/Users/jpaggi/Documents/combind/combind')
from containers import Protein
from shared_paths import proteins, shared_paths

def load_pdb(stats, helpers = 'pdb'):
    data = {}
    for protein in proteins:
        print(protein)
        fname = '{}/{}/scores/{}/summary/{}.tsv'.format(shared_paths['data'], protein, stats, helpers)
        with open(fname) as fp:
            fp.readline()
            for line in fp:
                print(line.strip().split('\t'))
                (mode, protein, ligand, alpha, features, _, combind_rmsd,
                 _, glide_rmsd, _, best_rmsd) = line.strip().split('\t')
                alpha, combind_rmsd,glide_rmsd, best_rmsd = float(alpha), float(combind_rmsd),float(glide_rmsd), float(best_rmsd)
                if mode not in data: data[mode] = {}
                if (alpha, features) not in data[mode]: data[mode][(alpha, features)] = {}
                if protein not in data[mode][(alpha, features)]: data[mode][(alpha, features)][protein] = {}
                data[mode][(alpha, features)][protein][ligand] = (float(combind_rmsd), float(glide_rmsd), float(best_rmsd))
        return data

def load_chembl(stats, helpers):
    data = {}
    for protein in proteins:
        print(protein)
        fname = '{}/{}/scores/{}/summary/{}.tsv'.format(shared_paths['data'], protein, stats, helpers)
        with open(fname) as fp:
            fp.readline()
            for line in fp:
                #if '\x00' in line: continue
                if 'None' in line: continue
                if len(line.split('\t')) != 12: print(line.split('\t'))
                (mode, protein, ligand, num_ligs, alpha, features, _, combind_rmsd,
                 _, glide_rmsd, _, best_rmsd) = line.strip().split('\t')
                alpha, combind_rmsd,glide_rmsd, best_rmsd = float(alpha), float(combind_rmsd),float(glide_rmsd), float(best_rmsd)
                num_ligs = int(num_ligs)
                if alpha == 0.5: continue
                if mode not in data: data[mode] = {}
                k = (num_ligs, alpha, features)
                if k not in data[mode]: data[mode][k] = {}
                if protein not in data[mode][k]: data[mode][k][protein] = {}
                data[mode][k][protein][ligand] = (float(combind_rmsd), float(glide_rmsd), float(best_rmsd))
    return data

def get_mcss_sizes(results):
    mcss_sizes = {}
    for protein_name, ligands in results.items():
        print(protein_name)
        protein = Protein(protein_name)
        lm = protein.lm
        lm.mcss.load_mcss()
        for ligand in ligands:
            crystal_lig = "{}_crystal_lig".format(lm.st)
            size = lm.mcss.get_mcss_size(ligand, crystal_lig)
            mcss_sizes[ligand] = size
    return mcss_sizes

##############################################################

class Marker:
    family = {
        'GPCR': ['5HT2B', 'A2AR', 'B1AR', 'B2AR', 'CHRM3','SMO', 'MGLUR5'],
        'Kinase': ['BRAF', 'CDK2', 'CHK1', 'JAK2', 'PLK1', 'MAPK14', 'MEK1'],
        'Transporter': ['SLC6A4', 'GLUT1', 'DAT', 'TRPV1'],
        'Nuclear Receptor': ['NR3C2', 'NR3C1', 'AR', 'VDR', 'ERA'],
        'Peptidase': ['F2', 'F10', 'F11', 'PLAU', 'P00760', 'BACE1'],
        'Other': ['PYGM', 'PTPN1', 'BRD4', 'HSP90AA1', 'PDE10A', 'SIGMAR1', 'ELANE']
    }

    markers = ['o', 'v', 's', '<', '>', 'X', 'D']
    colors  = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7"]

    def __init__(self):
        self.next = {fam: (self.colors[i], 0) for i, fam in enumerate(self.family)}

    def get_index(self, protein):
        fam = self.get_family(protein)
        return list(self.family.keys()).index(fam)

    def get_family(self, protein):
        for fam, prots in self.family.items():
            if protein in prots:
                return fam
        return 'Other'

    def __call__(self, protein):
        fam = self.get_family(protein)
        out = self.next[fam]
        self.next[fam] = (out[0], (out[1] + 1) % len(self.markers))
        return out[0], self.markers[out[1]]

def ligand_level_performance(results):
    x, y = [], []
    for prot, ligs in results.items():
        for lig, (combind, glide, best) in ligs.items():
            x += [glide]
            y += [combind]
    return x, y


def target_level_performance(results, valid_lig = lambda x, y: True, thresh = None):
    x, y, label = [], [], []
    for prot, ligs in results.items():
        _x, _y = [], []
        for lig, (combind, glide, best) in ligs.items():
            if valid_lig(prot, lig):
                if thresh is None:
                    _x += [glide]
                    _y += [combind]
                else:
                    _x += [glide <= thresh]
                    _y += [combind <= thresh]    
        if _x:
            label += [prot]
            x += [np.array(_x).mean()]
            y += [np.array(_y).mean()]
    return x, y, label

##################################################################

def ligand_level_plot(title, xlabel, ylabel, x, y, thresh):
    print(ttest_rel(x, y))
    print('{} improves pose for {} of {} ligands'.format(ylabel,
                                                         np.sum(np.array(x) > np.array(y)+.5),
                                                         len(x)))
    print('{} degrades pose for {} of {} ligands'.format(ylabel,
                                                         np.sum(np.array(x)+.5 < np.array(y)),
                                                         len(x)))

    f, ax = plt.subplots()
    plt.scatter(x, y, alpha = 0.5, s = 10)
    plt.xlabel(xlabel, fontsize = 16)
    plt.ylabel(ylabel, fontsize = 16)
    plt.title(title, fontsize = 20)
    plt.plot(range(int(math.ceil(max(x+y)))+1), linestyle='--', c = 'k')
    ax.set_aspect('equal', 'box')
    print('{} {}: {}, {}'.format(title, xlabel,
                                 sum(x) / float(len(x)),
                                 sum(np.array(x) <= thresh) /  float(len(x))))
    print('{} {}: {}, {}'.format(title, ylabel,
                                 sum(y) / float(len(y)),
                                 sum(np.array(y) <= thresh) /  float(len(y))))
    plt.show()

def target_level_plot(title, xlabel, ylabel, x, y, label):
    m = Marker()
    f, ax = plt.subplots()
    for i, (_x, _y, _label) in enumerate(zip(x, y, label)):
        color, marker = m(_label)
        plt.scatter(_x, _y, marker = marker, c = color, label = _label)
    plt.xlabel(xlabel, fontsize = 16)
    plt.ylabel(ylabel, fontsize = 16)
    plt.title(title, fontsize = 20)
    plt.plot(range(int(math.ceil(max(x+y)))+1), linestyle='--', c = 'k')
    ax.set_aspect('equal', 'box')
    print('{} {}: {}'.format(title, xlabel, sum(x) / float(len(x))))
    print('{} {}: {}'.format(title, ylabel, sum(y) / float(len(y))))
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles),
                                      key=lambda t: m.get_index(t[0])))
    plt.legend(handles, labels, bbox_to_anchor=(1.1, 1.05), ncol=3)
    plt.show()
    
######################################################################

def benchmark(results, thresh = 2.0, correct_only = False, families = None,
              mcss_sizes = None, low = None, high = None):
    if high is not None or low is not None:
        assert mcss_sizes is not None
    low  = -1 if low  is None else low
    high = 2  if high is None else high
    
    def valid_lig(prot, lig):
        def overlap():
            if low < 0 and high > 1:
                return True
            return low <= mcss_sizes[lig] < high
        
        def correct_exists():
            if not correct_only:
                return True
            return results[prot][lig][2] <= thresh
        
        def family():
            if families is None:
                return True
            return Marker().get_family(prot) in families

        return overlap() and correct_exists() and family()
    
    print('{} valid ligands'.format(sum(valid_lig(prot, lig)
                                        for prot, ligs in results.items()
                                        for lig in ligs)))

    # All ligands separately
    x, y = [], []
    for prot, ligs in results.items():
        for lig, (combind, glide, best) in ligs.items():
            if valid_lig(prot, lig):
                x += [glide]
                y += [combind]
    ligand_level_plot('All Ligands RMSD', 'Glide', 'ComBind', x, y, thresh)

    # By Protein (RMSD)
    x, y, label = target_level_performance(results, valid_lig, thresh = None)
    target_level_plot('Mean RMSD', 'Glide', 'ComBind', x, y, label)
    
    # By protein (% correct)
    x, y, label = target_level_performance(results, valid_lig, thresh = thresh)
    target_level_plot('Fraction Near-Native', 'Glide', 'ComBind', x, y, label)
