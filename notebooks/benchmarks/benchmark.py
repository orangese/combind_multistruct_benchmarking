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
        if protein == 'D2': continue
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
        'Transporter': ['SLC6A4', 'GLUT1', 'DAT'],
        'Nuclear Receptor': ['NR3C2', 'NR3C1', 'AR', 'VDR', 'ERA'],
        'Peptidase': ['F2', 'F10', 'F11', 'PLAU', 'P00760', 'BACE1'],
        'Other': ['PYGM', 'PTPN1', 'BRD4', 'HSP90AA1', 'PDE10A', 'SIGMAR1', 'ELANE', 'TRPV1', 'DHFR']
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

def ligand_level_performance(results, prots = [], correct_only = False):
    x, y, z = [], [], []
    for prot, ligs in results.items():
        if prots and prot not in prots: continue
        for lig, (combind, glide, best) in ligs.items():
            if best > 2.0: continue
            x += [glide]
            y += [combind]
            z += [best]
    return x, y, z


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
    f, ax = plt.subplots(figsize = (3, 10))
    
    
    order = ['GPCR', 'Transporter', 'Kinase', 'Nuclear Receptor', 'Peptidase', 'Other']
    order = order[::-1]
    
    def get_fam(prot):
        fam = [k for k, v in Marker.family.items() if prot in v]
        assert len(fam) == 1, prot
        return fam[0]
        
    
    def key(args):
        perf, _, prot = args
        fam = get_fam(prot)
        return order.index(fam)*1000 + perf
        
    
    x, y, label = zip(*sorted(zip(x, y, label), key = lambda args: key(args)))
    yticks = []
    i = 0
    last_fam = order[0]
    for _x, _y, prot in zip(x, y, label):
        if get_fam(prot) != last_fam:
            i += 2
            yticks += ['', '']
        yticks += [prot]
        plt.scatter(_y, i, c = 'g')
        plt.scatter(_x, i, c = 'm')
        last_fam = get_fam(prot)
        i += 1
    plt.yticks(range(len(yticks)), yticks, fontsize = 12)
    plt.title(title, fontsize = 20)
    print('{} {}: {}'.format(title, xlabel, sum(x) / float(len(x))))
    print('{} {}: {}'.format(title, ylabel, sum(y) / float(len(y))))
    plt.xlim(0, max(x+y))
    plt.show()
    
######################################################################

def benchmark(results, thresh = 2.0, correct_only = False, families = None,
              mcss_sizes = None, low = None, high = None, exclude = None):
    if high is not None or low is not None:
        assert mcss_sizes is not None
    low  = -1 if low  is None else low
    high = 2  if high is None else high
    
    def valid_lig(prot, lig):
        def not_excluded():
            if exclude is None:
                return True
            return prot not in exclude
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

        return overlap() and correct_exists() and family() and not_excluded()
    
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
