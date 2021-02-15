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

def load(input_csv, cwd, sim_cwd, xtal_cut=float('inf'), active_cut=float('inf')):
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
    df = df.loc[df['XTAL'] < xtal_cut]
    return df

################################################################################

# log-Adjusted ROC
# See https://pubmed.ncbi.nlm.nih.gov/20735049/
# Used recently in https://pubmed.ncbi.nlm.nih.gov/28760952/
def log_roc(x, y, lam=0.001):
    # Get rid of duplicated x values and those less than lambda.
    idx = np.array(list(x[:-1] != x[1:]) + [True])
    idx *= x >= lam
    x = x[idx]
    y = y[idx]

    # We need to have a value at lambda, so set a dummy one, if we don't.
    if x[0] != lam:
        y = np.array([lam*y[0]/x[0]] + list(y))
        x = np.array([lam]  + list(x))

    b = y[1:] - x[1:] * (y[1:] - y[:-1]) / (x[1:] - x[:-1])
    auc = np.sum((y[1:] - y[:-1]) / np.log(10))
    auc += np.sum(b * (np.log10(x[1:]) - np.log10(x[:-1])))
    auc /= np.log10(1/lam)
    auc -= .14462
    return auc, np.log10(x), y

def roc(x, y):
    return metrics.auc(x, y)

def ef(active, score, thresh=0.05):
    idx = np.argsort(-score)
    active = active[idx]

    thresh *= len(active)
    thresh = int(thresh)

    hit_rate = active[:thresh].mean()
    base_rate = active.mean() 
    return hit_rate / base_rate

############################################################################

def plot_cat(df, scores, colors, roc_ax=None, logroc_ax=None):
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

def plot_cat_all(xtal_cut, active_cut, helpers=5, helpers_sim=5, minimum=10,
                 weight_shape=10.0, weight_2D=50.0, weight_combind=1.0, plot=True,
                 scores = ['COMBIND', 'GLIDE', 'SHAPE_mean', '2D_mean', 'GLIDE+SHAPE', 'GLIDE+2D',
                           'COMBIND+2D', 'COMBIND+SHAPE'],
                 colors = {'COMBIND': 'g', 'GLIDE': 'b', 'SHAPE_mean':'m', '2D_mean': 'r',
                           'GLIDE+SHAPE': 'c', 'GLIDE+2D':'orange', 'COMBIND+2D': 'k',
                           'COMBIND+SHAPE': 'gray', 'COMBIND_GLIDE': 'lime', 'GLIDE_COMBIND': 'yellow'}):

    data = []
    for protein in os.listdir('/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind'):


        if plot:
            _, (roc_ax, logroc_ax) = plt.subplots(1, 2, figsize=(12, 5))
        else:
            roc_ax, logroc_ax = None, None
        aucs, logaucs, ef5s = [], [], []
        for i in range(5):
            cwd = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind/{}/scores/rd1_shape_{}/{}'.format(protein, helpers, i if helpers else 0)
            sim_cwd = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind/{}/scores/rd1_all_{}/{}'.format(protein, helpers_sim, i)
            input_csv = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_vs/DUDE/combind/{}/subset.smi'.format(protein)

            if not os.path.exists(cwd + '/combind.csv'): continue
            if not os.path.exists(cwd + '/glide.csv'): continue
            if not os.path.exists(cwd + '/shape.csv'): continue
            if not os.path.exists(cwd + '/similarity.csv'): continue
            if not os.path.exists(sim_cwd + '/similarity.csv'): continue

            
            stats_2D = '/oak/stanford/groups/rondror/users/jpaggi/2D_stats/'
            stats_SHAPE = '/oak/stanford/groups/rondror/users/jpaggi/SHAPE_stats/'
            nat_de = DensityEstimate.read('{}/nat_2D_{}_{}.de'.format(stats_2D, protein, helpers))
            ref_de = DensityEstimate.read('{}/ref_2D_{}_{}.de'.format(stats_2D, protein, helpers))
            nat_shape_de = DensityEstimate.read('{}/nat_SHAPE_{}_{}.de'.format(stats_SHAPE, protein, helpers))
            ref_shape_de = DensityEstimate.read('{}/ref_SHAPE_{}_{}.de'.format(stats_SHAPE, protein, helpers))

            df = load(input_csv, cwd, sim_cwd, xtal_cut=xtal_cut, active_cut=active_cut)
            
            df['COMBIND'] = (df['COMBIND']-df['COMBIND_GLIDE']) + weight_combind*df['COMBIND_GLIDE']
            
            df['GLIDE+SHAPE'] = df['GLIDE'] + weight_shape*df['SHAPE_mean']
            df['GLIDE+2D']    = df['GLIDE'] + weight_2D*df['2D_mean']
            
            df['COMBIND+SHAPE'] = df['COMBIND'] + weight_shape*df['SHAPE_mean']
            df['COMBIND+2D']    = df['COMBIND'] + weight_2D*df['2D_mean']

            df['Z:GLIDE+2D']    = Z(df['GLIDE']) + Z(df['2D_mean'])
            df['Z:COMBIND+2D']  = Z(df['COMBIND']) + Z(df['2D_mean'])

            df['Z:GLIDE+SHAPE']    = Z(df['GLIDE']) + Z(df['SHAPE_mean'])
            df['Z:COMBIND+SHAPE']  = Z(df['COMBIND']) + Z(df['SHAPE_mean'])

            df['Z:GLIDE+2D+SHAPE']    = Z(df['GLIDE']) + Z(df['SHAPE_mean']) + Z(df['2D_mean'])
            df['Z:COMBIND+2D+SHAPE']  = Z(df['COMBIND']) + Z(df['SHAPE_mean']) + Z(df['2D_mean'])
            df['Z:2D+SHAPE']  = Z(df['SHAPE_mean']) + Z(df['2D_mean'])

            df['N:2D_mean'] = np.log(nat_de(df['2D_mean'])) - np.log(ref_de(df['2D_mean']))
            df['N:SHAPE_mean'] = np.log(nat_shape_de(df['SHAPE_mean'])) - np.log(ref_shape_de(df['SHAPE_mean']))

            df['N:GLIDE+2D']    = df['GLIDE'] + df['N:2D_mean']
            df['N:COMBIND+2D']  = df['COMBIND'] + df['N:2D_mean']

            df['N:GLIDE+SHAPE']    = df['GLIDE'] + df['N:SHAPE_mean']
            df['N:COMBIND+SHAPE']  = df['COMBIND'] + df['N:SHAPE_mean']

            df['N:GLIDE+2D+SHAPE']    = df['GLIDE'] + df['N:SHAPE_mean'] + df['N:2D_mean']
            df['N:COMBIND+2D+SHAPE']  = df['COMBIND'] + df['N:SHAPE_mean'] + df['N:2D_mean']

            df['N:2D+SHAPE']  = df['N:SHAPE_mean'] + df['N:2D_mean']

            if sum(df['ACTIVE']) <= minimum:
                continue

            _aucs, _logaucs, _ef5s = plot_cat(df, scores, colors, roc_ax=roc_ax, logroc_ax=logroc_ax)
            aucs += [_aucs]
            logaucs += [_logaucs]
            ef5s += [_ef5s]

        if not aucs:
            if plot:
                plt.show()
            continue

        aucs = np.vstack(aucs).T
        logaucs = np.vstack(logaucs).T
        ef5s = np.vstack(ef5s).T

        if not np.any(np.isnan(logaucs.mean(axis=1))):
            data += [[protein] +  list(ef5s.mean(axis=1))]

        if plot:
            for score, auc, logauc in zip(scores, aucs, logaucs):
                label = '{} = {:.2f}'.format(score, np.mean(auc))
                roc_ax.scatter(-100, 0, s=50, marker='s', color=colors[score], label=label)
                label = '{} = {:.2f}'.format(score, np.mean(logauc))
                logroc_ax.scatter(-100, 0, s=50, marker='s', color=colors[score], label=label)

            roc_ax.scatter(-100, 0, s=0, label='# POS={}'.format(sum(df['ACTIVE'])))
            roc_ax.scatter(-100, 0, s=0, label='# NEG={}'.format(sum(~df['ACTIVE'])))
            logroc_ax.scatter(-100, 0, s=0, label='# POS={}'.format(sum(df['ACTIVE'])))
            logroc_ax.scatter(-100, 0, s=0, label='# NEG={}'.format(sum(~df['ACTIVE'])))


            roc_ax.set_ylim(0, 1)
            roc_ax.set_xlim(0, 1)
            roc_ax.legend(loc='lower right')
            roc_ax.set_ylabel('true positive rate', size=15)
            roc_ax.set_xlabel('false positive rate', size=15)
            logroc_ax.set_ylim(0, 1)
            logroc_ax.set_xlim(-3, 0)
            logroc_ax.legend(loc='upper left')
            logroc_ax.set_ylabel('true positive rate', size=15)
            logroc_ax.set_xlabel('log10 false positive rate', size=15)
            plt.suptitle(protein, size=20, y=0.96)
            plt.show()
    return pd.DataFrame(data, columns=['protein']+scores)

def aggregrate(data, metric1, metric2, xlabel, thresh=0.02):
    diffs = data[metric1] - data[metric2]

    win = np.sum(diffs > thresh)
    lose = np.sum(diffs < -thresh)
    tie = np.sum((diffs >= -thresh) & (diffs <= thresh))

    ax = plt.axes()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    bins = np.linspace(-0.5, 0.5, 36)
    plt.hist(diffs, bins=bins, edgecolor='b', color='white')
    plt.xlabel(xlabel.format(metric1, metric2), size=16)
    plt.ylabel('Count', size=16)
    plt.yticks(size=14)
    plt.xticks(size=14)
    plt.axvline(0, c='k', lw=3)

    y = ax.get_ylim()[1]
    ax.text(-.25, y, 'losses={}'.format(lose), fontsize=20,
            horizontalalignment='center', verticalalignment='bottom')

    ax.text(.25, y, 'wins={}'.format(win), fontsize=20,
            horizontalalignment='center', verticalalignment='bottom')
    plt.plot([-.5, -thresh], [y, y], c='k')
    plt.plot([thresh, .5], [y, y], c='k')
    plt.scatter(thresh, y, marker='|', s=100, c='k')
    plt.scatter(.5, y, marker='>', s=50, c='k')
    plt.scatter(-thresh, y, marker='|', s=100, c='k')
    plt.scatter(-.5, y, marker='<', s=50, c='k')
    plt.show()

