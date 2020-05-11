import click
import pandas as pd
import numpy as np
import os
from subprocess import run
from glob import glob
import scipy.stats
import sklearn.metrics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pandas as pd
import sys

@click.group()
def main():
    pass

@main.command()
@click.option('--n-train', default=10)
@click.option('--n-folds', default=7)
@click.option('--affinity-cut', default=1000)
@click.argument('input_csv')
@click.argument('root')
def setup(input_csv, root, n_train, n_folds, affinity_cut):
    np.random.seed(42)
    os.mkdir(root)
    df = pd.read_csv(input_csv)
    for i in range(n_folds):
        cwd = '{}/{}'.format(root, i)
        decoy_csv  = '{}/{}/decoy.csv'.format(root, i)
        binder_csv = '{}/{}/binder.csv'.format(root, i)
        train_csv  = '{}/{}/train.csv'.format(root, i)

        binders = df.loc[df['AFFINITY'] < affinity_cut]
        binders = binders.sample(n_train)

        decoys = df.loc[df['AFFINITY'] >= affinity_cut]
        decoys = decoys.sample(n_train)

        os.mkdir(cwd)
        decoys.to_csv(decoy_csv, index=False)
        binders.to_csv(binder_csv, index=False)
        pd.concat([decoys, binders]).to_csv(train_csv, index=False)

@main.command()
@click.option('--affinity-cut', default=1000)
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
def autoqsar(input_csv, root, affinity_cut):
    input_csv = os.path.abspath(input_csv)
    root = os.path.abspath(root)
    print(input_csv)
    for cwd in glob(root + '/[0-9]'):
        if os.path.exists('{}/autoqsar_preds.csv'.format(cwd)):
            continue
        run('$SCHRODINGER/utilities/autoqsar autoqsar.qzip -WAIT -build -i train.csv -y AFFINITY -cat -cuts {}'.format(affinity_cut),
            shell=True, cwd=cwd)
        run('$SCHRODINGER/utilities/autoqsar autoqsar.qzip -WAIT -test  -i {} -pred autoqsar_preds.csv'.format(input_csv),
            shell=True, cwd=cwd)

# for i in $(ls --color=none); do python ~/combind_dev/scripts/random_binders.py combind $i/subset.csv $i/scores/subset10_rd1 . $i; done;
@main.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
@click.argument('data')
@click.argument('protein')
def combind(input_csv, root, data, protein):
    input_csv = os.path.abspath(input_csv)
    root = os.path.abspath(root)
    data = os.path.abspath(data)
    for cwd in glob(root + '/[0-9]'):
        if not os.path.exists('{}/combind_preds.csv'.format(cwd)):
            run('$COMBINDHOME/main.py --data {} --ligands {}/binder.csv score {} all --pose-fname combind_poses.sc'.format(data, cwd, protein),
                shell=True, cwd=cwd)

        if not os.path.exists('{}/combind_preds.csv'.format(cwd)):
            run('$COMBINDHOME/main.py --data {} --ligands {} screen {} all --pose-fname combind_poses.sc --score-fname combind_preds.csv'.format(data, input_csv, protein),
                shell=True, cwd=cwd)

def get_fp(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

# for i in $(ls --color=none); do python ~/combind_dev/scripts/random_binders.py similarity $i/subset.csv $i/scores/subset10_rd1 . $i; done;
@main.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
@click.argument('data')
@click.argument('protein')
def similarity(input_csv, root, data, protein):
    ref = pd.read_csv('{}/{}/structures/pdb.csv'.format(data, protein))
    assert len(ref) == 1
    ref_smiles = ref.loc[0, 'SMILES']
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_fp = get_fp(ref_mol)

    
    df = pd.read_csv(input_csv)
    
    for cwd in glob(root + '/[0-9]'):
        active = pd.read_csv('{}/binder.csv'.format(cwd))
        active_fps = []
        for i, ligand in active.iterrows():
            mol = Chem.MolFromSmiles(ligand['SMILES'])
            active_fps += [get_fp(mol)]

        df['XTAL_sim'] = 0
        df['active_sim'] = 0
        for i, ligand in df.iterrows():
            try:
                mol = Chem.MolFromSmiles(ligand['SMILES'])
                fp = get_fp(mol)
                df.loc[i, 'XTAL_sim'] = DataStructs.TanimotoSimilarity(fp, ref_fp)
                df.loc[i, 'active_sim'] = max(DataStructs.TanimotoSimilarity(fp, active_fp)
                                              for active_fp in active_fps)
            except:
                print()
        df.to_csv('{}/similarity.csv'.format(cwd), index=False)
main()

def load_autoqsar(cwd):
    autoqsar = pd.read_csv(cwd+'/autoqsar_preds.csv')

    autoqsar['pred_cat'] = autoqsar['s_autoqsar_Pred_Class'] == '<1000.00'
    autoqsar['AUTOQSAR'] = 1 - autoqsar['r_autoqsar_Pred_Prob']
    autoqsar.loc[autoqsar['pred_cat'], 'AUTOQSAR'] = 1 - autoqsar.loc[autoqsar['pred_cat'], 'AUTOQSAR']
    autoqsar = autoqsar[['ID', 'AUTOQSAR']]
    autoqsar = autoqsar.set_index('ID')
    return autoqsar

def load_combind(cwd):
    combind = pd.read_csv(cwd+'/combind_preds.csv')
    combind['ID'] = [s.replace('_lig', '') for s in combind['ID']]
    combind = combind.set_index('ID')
    combind = - combind
    return combind

def load_affinity(input_csv):
    affinity = pd.read_csv(input_csv)
    affinity = affinity[['ID', 'CHEMBL', 'AFFINITY']]
    affinity = affinity.set_index('ID')
    return affinity

def load_train(cwd):
    train = pd.read_csv(cwd+'/train.csv')
    train = train['ID']
    return train

def load_similarity():
    similarity = pd.read_csv(cwd+'/similarity.csv')
    similarity = similarity[['ID', 'XTAL_sim', 'active_sim']].set_index('ID')
    return similarity

def load(cut=float('inf')):
    autoqsar = load_autoqsar()
    combind = load_combind()
    affinity = load_affinity()
    train = load_train()
    similarity = load_similarity()
    
    df = affinity.join(combind).join(autoqsar).join(similarity)
    df = df.fillna(0)
    df = df.loc[~df.index.isin(train)]
    df = df.loc[df['active_sim'] < cut]
    return df

def log_roc(ax, x, y, label='', lam=0.001):
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

    ax.plot(np.log10(x), y, label='{} = {:.2f}'.format(label, auc))
    
def roc(ax, x, y, label=''):
    auc = metrics.auc(x, y)
    ax.plot(fpr, tpr, label='{} = {:.2f}'.format(label, auc))

def plot(df, metrics=['GSCORE', 'CSCORE', 'AUTOQSAR']):
    print(df.shape)
    f, ax = plt.subplots(1, 2, figsize=(12, 5))

    for metric in metrics:
        fpr, tpr, _ = metrics.roc_curve(df['AFFINITY'] < 1000, df[metric])
        roc(ax[0], fpr, tpr, label=metric)
        log_roc(ax[1], fpr, tpr, label=metric)
    ax[0].legend()
    ax[1].legend()
    plt.show()

"""
# Build and then test qsar models using Schrodinger's autoqsar tool

# Categorical model
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -build -i train.csv -y AFFINITY -cat -cuts 1000
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -test  -i test.csv  -pred autoqsar_preds.csv

# Quantitative model
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -build -i train.csv -y AFFINITY -num -log -scale 1.0e-9
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -test  -i test.csv  -pred autoqsar_preds.csv
"""

def something():
    affinities = np.log10(affinities)-9
    gscores = np.array(gscores)
    cscores = np.array(cscores)
    xscores = np.array(xscores)
    ascores = np.array(ascores)

    ascores = np.minimum(xscores, ascores)

    queries = np.array(queries)[affinities < 0]
    gscores = gscores[affinities < 0]
    cscores = cscores[affinities < 0]
    xscores = xscores[affinities < 0]
    ascores = ascores[affinities < 0]
    affinities = affinities[affinities < 0]
    
    # evaluate performance for subsets of ligands with different
    # amounts of overlap with the XTAL ligand.
    for q in [0, 25, 50, 75]:
        mask = np.isin(np.array(queries), list(pose_cluster.keys()))
        mask |= np.array(xscores < np.quantile(xscores[affinities==-9], q/100))
        print(q, mask.sum())
        title = '{}_{}'.format(protein, q)
        evaluate_screen(title, affinities, gscores, cscores, xscores, ascores, mask)

    # and to closest active...
    for q in [0, 25, 50, 75]:
        mask = np.isin(np.array(queries), list(pose_cluster.keys()))
        mask |= np.array(ascores < np.quantile(ascores[affinities==-9], q/100))
        print(q, mask.sum())
        title = '{}_active_{}'.format(protein, q)
        evaluate_screen(title, affinities, gscores, cscores, xscores, ascores, mask)

    mask = np.isin(np.array(queries), list(pose_cluster.keys()))
    mask |= np.array(xscores < -0.4)
    print('xscore <= 0.4', mask.sum())
    title = '{}_{}'.format(protein, 4)
    evaluate_screen(title, affinities, gscores, cscores, xscores, ascores, mask)

    mask = np.isin(np.array(queries), list(pose_cluster.keys()))
    mask |= np.array(ascores < -0.4)
    print('ascore <= 0.4', mask.sum())
    title = '{}_active_{}'.format(protein, 4)
    evaluate_screen(title, affinities, gscores, cscores, xscores, ascores, mask)

def evaluate_screen(title, affinities, gscores, cscores, xscores, ascores, mask):
    def correlation(scores, plot_fname):
        rho, _ = scipy.stats.spearmanr(affinities[~mask], scores[~mask])
        r, _ = scipy.stats.pearsonr(affinities[~mask], scores[~mask])
        tau, _ = scipy.stats.kendalltau(affinities[~mask], scores[~mask])
        print(r, rho, tau)

        plt.scatter(affinities, scores, s=4, c='k')
        plt.scatter(affinities[mask], scores[mask], s=4, c='r')

        x = 0.25*plt.xlim()[1] + 0.75*plt.xlim()[0]
        y = 0.99*plt.ylim()[1] + 0.01*plt.ylim()[0]
        t = (r'$r={:.2f}$'
             "\n"
             r'$\rho={:.2f}$'
             "\n"
             r'$\tau={:.2f}$').format(r, rho, tau)
        plt.text(x, y, t, fontsize=20,
                 horizontalalignment='right', verticalalignment='top')
        plt.xlabel('log10 binding affinity', size=15)
        plt.ylabel('prediction', size=15)
        plt.title(plot_fname.split('.')[0])
        plt.savefig(plot_fname)
        plt.close()

    correlation(gscores, '{}_gscore.pdf'.format(title))
    correlation(cscores, '{}_cscore.pdf'.format(title))
    correlation(xscores, '{}_xscore.pdf'.format(title))
    correlation(ascores, '{}_ascore.pdf'.format(title))

    # ROC
    gfpr, gtpr, _ = sklearn.metrics.roc_curve(affinities[~mask]<-6, -gscores[~mask])
    cfpr, ctpr, _ = sklearn.metrics.roc_curve(affinities[~mask]<-6, -cscores[~mask])
    xfpr, xtpr, _ = sklearn.metrics.roc_curve(affinities[~mask]<-6, -xscores[~mask])
    afpr, atpr, _ = sklearn.metrics.roc_curve(affinities[~mask]<-6, -ascores[~mask])
    plt.plot(gfpr, gtpr, label='Glide   = {:.2f}'.format(sklearn.metrics.auc(gfpr, gtpr)))
    plt.plot(cfpr, ctpr, label='ComBind = {:.2f}'.format(sklearn.metrics.auc(cfpr, ctpr)))
    plt.plot(xfpr, xtpr, label='XTAL sim = {:.2f}'.format(sklearn.metrics.auc(xfpr, xtpr)))
    plt.plot(afpr, atpr, label='active sim = {:.2f}'.format(sklearn.metrics.auc(afpr, atpr)))

    plt.legend()
    plt.ylabel('true positive rate', size=15)
    plt.xlabel('false positive rate', size=15)
    plt.title(title)
    plt.savefig('{}_roc.pdf'.format(title))
    plt.close()

    # log-Adjusted ROC
    # See https://pubmed.ncbi.nlm.nih.gov/20735049/
    # Used recently in https://pubmed.ncbi.nlm.nih.gov/28760952/
    def log_roc(x, y, label='', lam=0.001):

        # Get rid of duplicated x values and those less than lambda.
        idx = np.array(list(x[:-1] != x[1:]) + [True])
        idx *= x >= lam
        x = x[idx]
        y = y[idx]

        # We need to have a value at lambda, so set a dummy one, if we don't.
        if x[0] != lam:
            x = np.array([lam]  + list(x))
            y = np.array([y[0]] + list(y))

        b = y[1:] - x[1:] * (y[1:] - y[:-1]) / (x[1:] - x[:-1])
        auc = np.sum((y[1:] - y[:-1]) / np.log(10))
        auc += np.sum(b * (np.log10(x[1:]) - np.log10(x[:-1])))
        auc /= np.log10(1/lam)

        plt.plot(np.log10(x), y, label='{} = {:.2f}'.format(label, auc))
    
    log_roc(gfpr, gtpr, 'Glide')
    log_roc(cfpr, ctpr, 'ComBind')
    log_roc(xfpr, xtpr, 'XTAL sim')
    log_roc(afpr, atpr, 'active sim')
    plt.legend()
    plt.ylabel('true positive rate', size=15)
    plt.xlabel('log10 false positive rate', size=15)
    plt.title(title)
    plt.savefig('{}_logroc.pdf'.format(title))
    plt.close()
