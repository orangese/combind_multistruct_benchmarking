import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
import scipy
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from score.density_estimate import DensityEstimate

def load_combind(cwd):
    df = pd.read_csv(cwd+'/combind.csv')
    df = df.set_index('ID')
    df = df.rename(columns={'GLIDE': 'COMBIND_GLIDE'})
    df['COMBIND_GLIDE'] *= -1
    df[df < 0] = 0
    return df

def load_glide(cwd):
    df = pd.read_csv(cwd+'/glide.csv')
    df = df.set_index('ID')
    df = df.rename(columns={'COMBIND': 'GLIDE_COMBIND'})
    df['GLIDE'] *= -1
    df[df < 0] = 0
    return df

def load_shape(cwd):
    df = pd.read_csv(cwd+'/shape.csv')
    df = df.set_index('ID')
    return df

def load_similarity(cwd):
    df = pd.read_csv(cwd+'/similarity.csv')
    df = df.set_index('ID')
    df = df.rename(columns={'XTAL_sim': 'XTAL', 'active_sim_mean': '2D_mean', 'active_sim_max': '2D_max'})
    df = df[['XTAL', '2D_mean', '2D_max']]
    return df

def load(input_csv, cwd, sim_cwd, xtal_cut=float('inf'), active_cut=float('inf'), active_min_cut=0):
    df = pd.read_csv(input_csv, sep=' ')
    df['ACTIVE'] = df['ID'].str.contains('CHEMBL')
    df = df.set_index('ID')

    combind = load_combind(cwd)
    glide = load_glide(cwd)
    similarity = load_similarity(cwd)
    shape = load_shape(cwd)
    
    similarity_cutoff = load_similarity(sim_cwd)
    similarity_cutoff['sim_2D_max'] = similarity_cutoff['2D_max']
    similarity_cutoff = similarity_cutoff[['sim_2D_max']]

    df = df.join(combind).join(glide).join(similarity).join(shape).join(similarity_cutoff)

    df = df.fillna(0)
    df = df.loc[df['sim_2D_max'] < active_cut]
    df = df.loc[df['sim_2D_max'] >= active_min_cut]
    df = df.loc[df['XTAL'] < xtal_cut]
    return df

################################################################################

# log-Adjusted ROC
# See https://pubmed.ncbi.nlm.nih.gov/20735049/
# Used recently in https://pubmed.ncbi.nlm.nih.gov/28760952/
def log_roc(x, y, lam=0.001):
    # We need to have a value at lambda.
    for i in range(len(x)):
        if x[i] > lam:
            m = (y[i]-y[i-1]) / (x[i]-x[i-1])
            b = y[i] - m*x[i]
            x = np.array([lam] + list(x[i:]))
            y = np.array([m*lam + b] + list(y[i:]))
            break
        elif x[i] == lam:
            x = x[i:]
            y = y[i:]
            break
    else:
        assert False
    
    # Get rid of duplicated x values.
    idx = np.array(list(x[:-1] != x[1:]) + [True])
    x = x[idx]
    y = y[idx]

    b = y[1:] - x[1:] * (y[1:] - y[:-1]) / (x[1:] - x[:-1])
    auc = np.sum((y[1:] - y[:-1]) / np.log(10))
    auc += np.sum(b * (np.log10(x[1:]) - np.log10(x[:-1])))
    auc /= np.log10(1/lam)
    auc -= .14462
    return auc, np.log10(x), y

def roc(x, y):
    return metrics.auc(x, y)

def ef(active, score, thresh=0.01):
    idx = np.argsort(-score)
    active = active[idx]

    thresh *= len(active)
    thresh = int(thresh)

    hit_rate = active[:thresh].mean()
    base_rate = active.mean()
    return hit_rate / base_rate

############################################################################

def plot_cat(df, scores, colors=None, roc_ax=None, logroc_ax=None):
    aucs, logaucs, ef5s = [], [], []
    for score in scores:
        fpr, tpr, _ = metrics.roc_curve(df['ACTIVE'], df[score])
        auc = roc(fpr, tpr)
        log_auc, log_fpr, log_tpr = log_roc(fpr, tpr)
        ef5 = ef(df['ACTIVE'].to_numpy(), df[score].to_numpy())
        
        aucs += [auc]
        logaucs += [log_auc]
        ef5s += [ef5]
        
        if roc_ax:
            label = '{} = {:.2f}'.format(score, auc)
            color = colors[score]
            roc_ax.plot(fpr, tpr, color=color)
           
        if logroc_ax:
            label = '{} = {:.2f}'.format(score, log_auc)
            color = colors[score]
            logroc_ax.plot(log_fpr, log_tpr, color=color)
    return aucs, logaucs, ef5s

def Z(x):
    return (x - x.mean()) / np.sqrt(x.var())

def R(x):
    order = np.argsort(-x)
    ranks = np.argsort(order)+1
    return -ranks

def plot_cat_all(xtal_cut, active_cut, helpers=5, helpers_sim=5,
                 active_min_cut=0, minimum_active=10, minimum_decoy=100,
                 metric='logroc', scores=None, minimum_repeats=1,
                 stats_2D='/oak/stanford/groups/rondror/users/jpaggi/combindvs_initial_submission/2D_stats/',
                 stats_SHAPE='/oak/stanford/groups/rondror/users/jpaggi/combindvs_initial_submission/SHAPE_stats/',
                 root='/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind'):

    data = []
    for protein in os.listdir(root):
        print(protein)
        aucs, logaucs, ef5s = [], [], []
        for i in range(5):
            cwd = '{}/{}/scores/all_rd1_shape_{}/{}'.format(root, protein, helpers, i if helpers else 0)
            sim_cwd = '{}/{}/scores/all_rd1_shape_{}/{}'.format(root, protein, helpers_sim, i)
            input_csv = '{}/{}/all.smi'.format(root, protein)

            if not os.path.exists(cwd + '/combind.csv'): continue
            if not os.path.exists(cwd + '/glide.csv'): continue
            if not os.path.exists(cwd + '/shape.csv'): continue
            if not os.path.exists(cwd + '/similarity.csv'): continue
            if not os.path.exists(sim_cwd + '/similarity.csv'): continue
            
            nat_de = DensityEstimate.read('{}/nat_2D_{}_{}.de'.format(stats_2D, protein, helpers))
            ref_de = DensityEstimate.read('{}/ref_2D_{}_{}.de'.format(stats_2D, protein, helpers))
            nat_shape_de = DensityEstimate.read('{}/nat_SHAPE_{}_{}.de'.format(stats_SHAPE, protein, helpers))
            ref_shape_de = DensityEstimate.read('{}/ref_SHAPE_{}_{}.de'.format(stats_SHAPE, protein, helpers))

            df = load(input_csv, cwd, sim_cwd, xtal_cut=xtal_cut,
                      active_cut=active_cut, active_min_cut=active_min_cut)

            if sum(df['ACTIVE']) <= minimum_active:
                continue

            if len(df)-sum(df['ACTIVE']) <= minimum_decoy:
                continue

            df['BEST'] = df['ACTIVE'].astype(int)

            # Z-score
            df['Z:COMBIND'] = Z(df['COMBIND'])
            df['Z:GLIDE'] = Z(df['GLIDE'])
            df['Z:2D'] = Z(df['2D_mean'])
            df['Z:SHAPE'] = Z(df['SHAPE_mean'])
            
            df['Z:GLIDE+2D'] = df['Z:GLIDE'] + df['Z:2D']
            df['Z:COMBIND+2D'] = df['Z:COMBIND'] + df['Z:2D']
            df['Z:GLIDE+SHAPE'] = df['Z:GLIDE'] + df['Z:SHAPE']
            df['Z:COMBIND+SHAPE'] = df['Z:COMBIND'] + df['Z:SHAPE']
            df['Z:GLIDE+2D+SHAPE'] = df['Z:GLIDE'] + df['Z:SHAPE'] + df['Z:2D']
            df['Z:COMBIND+2D+SHAPE'] = df['Z:COMBIND'] + df['Z:SHAPE'] + df['Z:2D']
            df['Z:2D+SHAPE'] = df['Z:SHAPE'] + df['Z:2D']

            # Product rank.
            df['R:COMBIND'] = R(df['COMBIND'])
            df['R:GLIDE'] = R(df['GLIDE'])
            df['R:2D'] = R(df['2D_mean'])
            df['R:SHAPE'] = R(df['SHAPE_mean'])

            df['R:GLIDE+2D'] = -df['R:GLIDE'] * df['R:2D']
            df['R:COMBIND+2D'] = -df['R:COMBIND'] * df['R:2D']
            df['R:GLIDE+SHAPE'] = -df['R:GLIDE'] * df['R:SHAPE']
            df['R:COMBIND+SHAPE'] = -df['R:COMBIND'] * df['R:SHAPE']
            df['R:GLIDE+2D+SHAPE'] = df['R:GLIDE'] * df['R:2D'] * df['R:SHAPE']
            df['R:COMBIND+2D+SHAPE'] = df['R:COMBIND'] * df['R:2D'] * df['R:SHAPE']
            df['R:2D+SHAPE'] =  -df['R:2D'] * df['R:SHAPE']

            # Naive bayes.
            df['N:2D'] = np.log(nat_de(df['2D_mean'])) - np.log(ref_de(df['2D_mean']))
            df['N:SHAPE'] = np.log(nat_shape_de(df['SHAPE_mean'])) - np.log(ref_shape_de(df['SHAPE_mean']))

            df['N:GLIDE+2D'] = df['GLIDE'] + df['N:2D']
            df['N:COMBIND+2D'] = df['COMBIND'] + df['N:2D']

            df['N:GLIDE+SHAPE'] = df['GLIDE'] + df['N:SHAPE']
            df['N:COMBIND+SHAPE'] = df['COMBIND'] + df['N:SHAPE']

            df['N:GLIDE+2D+SHAPE'] = df['GLIDE'] + df['N:SHAPE'] + df['N:2D']
            df['N:COMBIND+2D+SHAPE'] = df['COMBIND'] + df['N:SHAPE'] + df['N:2D']

            df['N:2D+SHAPE'] = df['N:SHAPE'] + df['N:2D']

            _aucs, _logaucs, _ef5s = plot_cat(df, scores)
            aucs += [_aucs]
            logaucs += [_logaucs]
            ef5s += [_ef5s]

        if len(logaucs) < minimum_repeats:
            continue

        aucs = np.vstack(aucs).T
        logaucs = np.vstack(logaucs).T
        ef5s = np.vstack(ef5s).T

        if metric == 'logroc':
            _data = logaucs
        elif metric == 'ef1':
            _data = ef5s
        elif metric == 'roc':
            _data = aucs
        else:
            assert False, metric
        data += [[protein] + list(_data.mean(axis=1))]
    return pd.DataFrame(data, columns=['protein']+scores)
