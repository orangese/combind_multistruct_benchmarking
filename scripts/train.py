import numpy as np
import sys
import pandas as pd
import click
from glob import glob
import pickle
from pygam import LogisticGAM, s, l, f
from sklearn.linear_model import LogisticRegression
from density_estimate import DensityEstimate

@click.group()
def main():
    pass

@main.command()
@click.argument('input_csv')
@click.argument('out_root')
@click.argument('sd', default=0.04)
@click.argument('weighting', default='unweighted')
def train_nb(input_csv, out_root, sd, weighting):
    data = pd.read_csv(input_csv)
    data = data.set_index(['protein', 'ligand1', 'ligand2'])

    if weighting == 'unweighted':
        data['W'] = 1
    elif weighting == 'ligand_pair':
        data['W'] = data['W_lig_pair']
    elif weighting == 'ligand':
        data['W'] = data['W_lig_pair']*data['W_ligand1']
    elif weighting == 'protein':
        data['W'] = data['W_lig_pair']*data['W_ligand1']*data['W_ligand2']

    for feature in ['hbond', 'saltbridge', 'contact', 'mcss', 'shape']:
        if feature == 'mcss':
            domain = (-6, 0)
            _sd = sd*6
            _data = data.loc[data.no_mcss == 0]
        else:
            domain = (0, 1)
            _sd = sd
            _data = data

        print(_data['W'].to_numpy().shape)

        nat = DensityEstimate(domain=domain, sd=sd, reflect=False).fit(_data.loc[_data['native'], feature], _data.loc[_data['native'], 'W'].to_numpy())
        nat.write('{}/nb_{}_nat.de'.format(out_root, feature))

        ref = DensityEstimate(domain=domain, sd=sd, reflect=False).fit(_data.loc[~_data['native'], feature], _data.loc[~_data['native'], 'W'].to_numpy())
        ref.write('{}/nb_{}_ref.de'.format(out_root, feature))

@main.command()
@click.argument('input_csv')
@click.argument('out_root')
@click.argument('features', default='gscore,shape,mcss,hbond,saltbridge,contact,no_mcss')
@click.argument('weighting', default='unweighted')
def train_lr(input_csv, out_root, features, weighting):
    features = features.split(',')
    data = pd.read_csv(input_csv)
    data = data.set_index(['protein', 'ligand1', 'ligand2'])

    if weighting == 'unweighted':
        data['W'] = 1
    elif weighting == 'ligand_pair':
        data['W'] = data['W_lig_pair']
    elif weighting == 'ligand':
        data['W'] = data['W_lig_pair']*data['W_ligand1']
    elif weighting == 'protein':
        data['W'] = data['W_lig_pair']*data['W_ligand1']*data['W_ligand2']

    lr = LogisticRegression(solver='lbfgs').fit(data[features], data['native'], data['W'])

    pickle.dump(lr, open('{}/lr_{}.p'.format(out_root, '-'.join(features)), "wb"))

@main.command()
@click.argument('input_csv')
@click.argument('out_root')
@click.argument('features', default='gscore,shape,mcss,hbond,saltbridge,contact,no_mcss')
@click.argument('weighting', default='unweighted')
def train_lg(input_csv, out_root, features, weighting):
    features = features.split(',')
    data = pd.read_csv(input_csv)
    data = data.set_index(['protein', 'ligand1', 'ligand2'])

    if weighting == 'unweighted':
        data['W'] = 1
    elif weighting == 'ligand_pair':
        data['W'] = data['W_lig_pair']
    elif weighting == 'ligand':
        data['W'] = data['W_lig_pair']*data['W_ligand1']
    elif weighting == 'protein':
        data['W'] = data['W_lig_pair']*data['W_ligand1']*data['W_ligand2']

    data['W'] /= data['W'].mean()

    for i, feature in enumerate(features):
        if feature == 'gscore':
            if i:
                terms += l(i)
            else:
                terms = l(i)
        elif 'no_' in feature:
            if i:
                terms += f(i)
            else:
                terms = f(i)
        else:
            if i:
                terms += s(i)
            else:
                terms = s(i)

    lg = LogisticGAM(terms)
    lg.fit(data[features], data['native'], weights=data['W'].to_numpy())

    pickle.dump(lg, open('{}/lg_{}.p'.format(out_root, '-'.join(features)), "wb"))

main()