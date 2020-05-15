import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
import scipy
import numpy as np
import os

def load_autoqsar_cat(cwd, positive='<1000.00'):
    autoqsar = pd.read_csv(cwd+'/autoqsar_preds.csv')

    autoqsar['pred_cat'] = autoqsar['s_autoqsar_Pred_Class'] == positive
    autoqsar['AUTOQSAR'] = autoqsar['r_autoqsar_Pred_Prob']
    autoqsar.loc[~autoqsar['pred_cat'], 'AUTOQSAR'] = \
                        1 - autoqsar.loc[~autoqsar['pred_cat'], 'AUTOQSAR']
    autoqsar = autoqsar[['ID', 'AUTOQSAR']]
    autoqsar = autoqsar.set_index('ID')
    return autoqsar

def load_autoqsar_num(cwd):
    autoqsar = pd.read_csv(cwd+'/autoqsar_preds.csv')
    autoqsar['AUTOQSAR'] = - autoqsar['r_autoqsar_Pred_Y']
    autoqsar = autoqsar[['ID', 'AUTOQSAR']]
    autoqsar = autoqsar.set_index('ID')
    return autoqsar

def load_combind(cwd):
    combind = pd.read_csv(cwd+'/combind_preds.csv')
    combind['ID'] = [s.replace('_lig', '') for s in combind['ID']]
    combind = combind.set_index('ID')
    combind.loc[:, 'CSCORE'] *= -1
    combind.loc[:, 'CSCORE_GPOSE'] *= -1
    combind[combind > 0] = 0
    return combind

def load_affinity(input_csv):
    affinity = pd.read_csv(input_csv)
    cols = ['ID', 'CHEMBL', 'AFFINITY', 'LOGAFFINITY']
    cols = [col for col in cols if col in affinity]
    affinity = affinity[cols]
    affinity = affinity.set_index('ID')
    return affinity

def load_train(cwd, train_fname):
    train = pd.read_csv(cwd+'/'+train_fname)
    train = train['ID']
    return train

def load_similarity(cwd):
    similarity = pd.read_csv(cwd+'/similarity.csv')
    if 'decoy_sim' not in similarity:
        similarity['decoy_sim'] = 0
    similarity = similarity[['ID', 'XTAL_sim', 'active_sim', 'decoy_sim']].set_index('ID')
    return similarity

def load(input_csv, cwd, xtal_cut=float('inf'), active_cut=float('inf'),
         decoy_cut=float('inf'), train_fname='train.csv', cat=None, num=False):
    if num:
        assert cat is None
        autoqsar = load_autoqsar_num(cwd)
    else:
        assert cat is not None
        autoqsar = load_autoqsar_cat(cwd, cat)

    combind = load_combind(cwd)
    affinity = load_affinity(input_csv)
    train = load_train(cwd, train_fname)
    similarity = load_similarity(cwd)

    df = affinity.join(combind).join(autoqsar).join(similarity)
    df = df.fillna(0)
    df['train'] = df.index.isin(train)
    df = df.loc[df['active_sim'] < active_cut]
    df = df.loc[df['XTAL_sim'] < xtal_cut]
    df = df.loc[df['decoy_sim'] < decoy_cut]
    return df

################################################################################

# log-Adjusted ROC
# See https://pubmed.ncbi.nlm.nih.gov/20735049/
# Used recently in https://pubmed.ncbi.nlm.nih.gov/28760952/
def log_roc(ax, x, y, label=None, lam=0.001, color=None):
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

    if label != None:
        label = '{} = {:.2f}'.format(label, auc)

    ax.plot(np.log10(x), y, label=label, color=color)
    return auc
    
def roc(ax, x, y, label=None, color=None):
    auc = metrics.auc(x, y)

    ax.plot(x, y, label=label, color=color)
    return auc

def plot_cat(df, ax, scores, colors):
    aucs, logaucs = [], []
    for score in scores:
        fpr, tpr, _ = metrics.roc_curve(df['AFFINITY'] < 1000, df[score])
        aucs += [roc(ax[0], fpr, tpr, color=colors[score])]
        logaucs += [log_roc(ax[1], fpr, tpr, color=colors[score])]
    return aucs, logaucs

def plot_cat_all(xtal_cut, active_cut, decoy_cut,
                 scores = ['CSCORE', 'GSCORE', 'AUTOQSAR'],
                 colors = {'GSCORE': 'b', 'CSCORE': 'g', 'AUTOQSAR':'m'},
                 names = {'GSCORE': 'Glide', 'CSCORE': 'ComBind', 'AUTOQSAR':'AUTOQSAR'}):
    data = []
    for protein in os.listdir('/Users/jpaggi/sherlock/oak/users/jpaggi/VS/DUDE/combind'):
        f, ax = plt.subplots(1, 2, figsize=(12, 5))
        aucs, logaucs = [], []
        for i in range(7):
            cwd = '/Users/jpaggi/sherlock/oak/users/jpaggi/VS/DUDE/combind/{}/scores/subset10_rd1/{}'.format(protein, i)
            input_csv = '/Users/jpaggi/sherlock/oak/users/jpaggi/VS/DUDE/combind/{}/subset.csv'.format(protein)

            if not os.path.exists(cwd + '/combind_preds.csv'): continue
            if not os.path.exists(cwd + '/autoqsar_preds.csv'): continue

            df = load(input_csv, cwd, cat='<1000.00', xtal_cut=xtal_cut, active_cut=active_cut, decoy_cut=decoy_cut)
            df.loc[:, 'GSCORE'] *= -1
            df.loc[:, 'CSCORE'] *= -1
            df.loc[:, 'GSCORE_CPOSE'] *= -1
            df.loc[:, 'CSCORE_GPOSE'] *= -1
            _aucs, _logaucs = plot_cat(df.loc[~df['train']], ax, scores, colors)
            aucs += [_aucs]
            logaucs += [_logaucs]

        if not aucs:
            plt.show()
            continue

        aucs = np.vstack(aucs).T
        logaucs = np.vstack(logaucs).T

        if not np.any(np.isnan(logaucs.mean(axis=1))):
            data += [[protein] + list(logaucs.mean(axis=1))]

        for score, auc, logauc in zip(scores, aucs, logaucs):
            label = '{} = {:.2f}'.format(names[score], np.mean(auc))
            ax[0].scatter(-100, 0, s=50, marker='s', color=colors[score], label=label)
            label = '{} = {:.2f}'.format(names[score], np.mean(logauc))
            ax[1].scatter(-100, 0, s=50, marker='s', color=colors[score], label=label)

        for i in [0, 1]:
            ax[i].scatter(-100, 0, s=0, label='# POS={}'.format(sum(df['AFFINITY'] < 1000)))
            ax[i].scatter(-100, 0, s=0, label='# NEG={}'.format(sum(df['AFFINITY'] >= 1000)))

        ax[0].set_ylim(0, 1)
        ax[0].set_xlim(0, 1)
        ax[0].legend(loc='lower right')
        ax[0].set_ylabel('true positive rate', size=15)
        ax[0].set_xlabel('false positive rate', size=15)
        ax[1].set_ylim(0, 1)
        ax[1].set_xlim(-3, 0)
        ax[1].legend(loc='upper left')
        ax[1].set_ylabel('true positive rate', size=15)
        ax[1].set_xlabel('log10 false positive rate', size=15)
        plt.suptitle(protein, size=20, y=0.96)
        plt.show()
    return pd.DataFrame(data, columns=['protein']+[names[s] for s in scores])

################################################################################

def plot_num(ax, mask, affinities, scores, plot=True):
    if plot:
        ax.scatter(affinities, scores, s=4, c='k')
        ax.scatter(affinities[mask], scores[mask], s=4, c='r')
    
    tau, tau_p = scipy.stats.kendalltau(affinities[~mask], scores[~mask])
    r, r_p = scipy.stats.pearsonr(affinities[~mask], scores[~mask])
    
    if plot:
        x = 0.35*ax.get_xlim()[1] + 0.65*ax.get_xlim()[0]
        y = 0.99*ax.get_ylim()[1] + 0.01*ax.get_ylim()[0]
        t = (r'$\tau={:.2f}$'
             "\n"
             r'$r={:.2f}$').format(tau, r)
        ax.text(x, y, t, fontsize=20,
               horizontalalignment='right', verticalalignment='top')
        ax.set_xlabel('log10 binding affinity', size=15)
    return tau, r

def plot_num_all(affinity_cut, xtal_cut, active_cut, plot=True):
    data = []
    for protein in os.listdir('/Users/jpaggi/sherlock/oak/users/jpaggi/VS/CHEMBL'):
        taus, rs = [], []
        for i in range(7):
            cwd = '/Users/jpaggi/sherlock/oak/users/jpaggi/VS/CHEMBL/{}/scores/rd1/{}'.format(protein, i)
            if not os.path.exists(cwd + '/combind_preds.csv'): continue
            if not os.path.exists(cwd + '/autoqsar_preds.csv'): continue

            df = load('/Users/jpaggi/sherlock/oak/users/jpaggi/VS/CHEMBL/{}/chembl/chosen.csv'.format(protein),
                      cwd, xtal_cut=xtal_cut, active_cut=active_cut, train_fname='binder.csv',
                      num=True)
            df = df.loc[df['LOGAFFINITY'] < affinity_cut]
            
            if sum(~df['train']) < 30:
                continue

            if plot:
                f, ax = plt.subplots(1, 3, figsize=(16, 4))
                ax[0].set_ylabel('CSCORE', size=15)
                ax[1].set_ylabel('GSCORE', size=15)
                ax[2].set_ylabel('AUTOQSAR', size=15)
                ax[1].set_title('{}_{}'.format(protein, i), size=18)
            else:
                ax = [0, 0, 0]

            ctau, cr = plot_num(ax[0], df['train'], df['LOGAFFINITY'], df['CSCORE'], plot)
            gtau, gr = plot_num(ax[1], df['train'], df['LOGAFFINITY'], df['GSCORE'], plot)
            atau, ar = plot_num(ax[2], df['train'], df['LOGAFFINITY'], df['AUTOQSAR'], plot)
            taus += [[ctau, gtau, atau]]
            rs += [[cr, gr, ar]]

            if plot:
                plt.show()

        if not taus:
            continue

        taus = np.vstack(taus)
        rs = np.vstack(taus)
        data += [[protein] + list(taus.mean(axis=0))]
    return pd.DataFrame(data, columns=['protein', 'ComBind', 'Glide', 'AUTOQSAR'])

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
