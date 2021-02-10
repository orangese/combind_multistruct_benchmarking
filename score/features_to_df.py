import numpy as np
import sys
import pandas as pd
from features.features import Features
import click
from glob import glob
from utils import mp

def load_protein(data, protein):
    interactions = ['hbond',  'saltbridge', 'contact', 'shape', 'mcss']
    features = Features(data + '/' + protein, max_poses=100)
    features.load_features(interactions)
    features = features.raw

    # get cross-docked ligands
    ligands = []
    for ligand in sorted(features['gscore'].keys()):
        if 'native' in ligand: continue
        lig, grid = ligand.split('-to-')
        if '_lig' in lig:
            lig = lig.replace('_lig', '')
        if lig != grid and 'CHEMBL' not in lig:
            ligands += [ligand]

    df = []
    for i, ligand1 in enumerate(ligands):
        for ligand2 in ligands[i+1:]:
            for r1 in range(len(features['gscore'][ligand1])):
                for r2 in range(len(features['gscore'][ligand2])):
                    feats = [features[interaction][(ligand1, ligand2)][r1, r2]
                             for interaction in interactions]
                    gscore1 = features['gscore'][ligand1][r1]
                    gscore2 = features['gscore'][ligand2][r2]
                    rmsd1 = features['rmsd'][ligand1][r1]
                    rmsd2 = features['rmsd'][ligand2][r2]
                    df += [[protein,
                            ligand1, ligand2,
                            r1, r2,
                            gscore1, gscore2,
                            rmsd1, rmsd2]
                            + feats]
    return pd.DataFrame(df, columns=['protein',
                                     'ligand1', 'ligand2',
                                     'rank1', 'rank2',
                                     'gscore1', 'gscore2',
                                     'rmsd1', 'rmsd2']
                                     +interactions)

def load_protein_top(data, protein, n_native, n_decoy):
    interactions = ['hbond',  'saltbridge', 'contact', 'mcss', 'shape', 'pipi', 'pi-t']
    features = Features(data + '/' + protein, max_poses=100)
    features.load_features(interactions)
    features = features.raw

    # get cross-docked ligands
    ligands = []
    for ligand in sorted(features['gscore'].keys()):
        if 'native' in ligand: continue
        lig, grid = ligand[:-3].split('-to-')
        if lig != grid:
            ligands += [ligand]

    df = []
    for ligand1 in ligands:
        for r1 in range(len(features['gscore'][ligand1])):
            if features['rmsd'][ligand1][r1] <= 2.05:
                break
        else:
            continue

        for ligand2 in ligands:
            if ligand2 == ligand1: continue
    
            r2_native = [r for r in range(len(features['gscore'][ligand2]))
                         if features['rmsd'][ligand2][r] <= 2.05]
            r2_decoy = [r for r in range(len(features['gscore'][ligand2]))
                         if features['rmsd'][ligand2][r] > 2.05]
            r2_native = r2_native[:n_native]
            r2_decoy = r2_decoy[:n_decoy]
            for r2 in r2_native + r2_decoy:

                if ligand1 < ligand2:
                    feats = [features[interaction][(ligand1, ligand2)][r1, r2]
                             for interaction in interactions]
                else:
                    feats = [features[interaction][(ligand2, ligand1)][r2, r1]
                             for interaction in interactions]
                
                gscore1 = features['gscore'][ligand1][r1]
                gscore2 = features['gscore'][ligand2][r2]
                rmsd1 = features['rmsd'][ligand1][r1]
                rmsd2 = features['rmsd'][ligand2][r2]
                df += [[protein,
                        ligand1, ligand2,
                        r1, r2,
                        gscore1, gscore2,
                        rmsd1, rmsd2]
                        + feats]
    return pd.DataFrame(df, columns=['protein',
                                     'ligand1', 'ligand2',
                                     'rank1', 'rank2',
                                     'gscore1', 'gscore2',
                                     'rmsd1', 'rmsd2']
                                     +interactions)

@click.group()
def main():
    pass

@main.command()
@click.argument('protein')
@click.argument('data', default='/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020_v2')
@click.argument('output_root', default='pairs')
def run(protein, data, output_root):
    df = load_protein(data, protein)
    df.to_csv('{}/{}.csv'.format(output_root, protein), index=False)

@main.command()
@click.argument('protein')
@click.argument('data', default='/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020_v2')
@click.argument('output_root', default='pairs')
@click.option('--n-native', default=1)
@click.option('--n-decoy', default=1)
def run_top(protein, data, output_root,  n_native, n_decoy):
    df = load_protein_top(data, protein,  n_native, n_decoy)
    df.to_csv('{}/{}.csv'.format(output_root, protein), index=False)

@main.command()
@click.argument('input_pattern')
@click.argument('output_csv')
def merge(input_pattern, output_csv):
    fnames = glob(input_pattern)
    df = pd.concat([pd.read_csv(fname) for fname in fnames])
    df.to_csv(output_csv, index=False)

@main.command()
@click.argument('input_csv')
@click.argument('output_csv')
def transform(input_csv, output_csv):
    data = pd.read_csv(input_csv)
    data = data.set_index(['protein', 'ligand1', 'ligand2'])
    data['gscore1'] = -data['gscore1']
    data['gscore2'] = -data['gscore2']
    data['mcss']    = -data['mcss']
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')

    data['no_mcss'] = data.mcss == -float('inf')
    data.loc[data.no_mcss, 'mcss'] = 0.0
    data.loc[(data['mcss'] < -6).to_numpy(), 'mcss'] = -6

    data['no_saltbridge'] = 0
    data.loc[data.groupby(level=[0, 1, 2]).max()['saltbridge'] == 0.5, 'no_saltbridge'] = 1

    data_rev = data.reset_index().rename(columns={'gscore1': 'gscore2', 'gscore2': 'gscore1',
                                      'rank1': 'rank2', 'rank2': 'rank1',
                                      'rmsd1': 'rmsd2', 'rmsd2': 'rmsd1',
                                      'ligand1': 'ligand2', 'ligand2': 'ligand1'})
    data_rev = data_rev.set_index(['protein', 'ligand1', 'ligand2']).sort_index()
    data_rev = data_rev[list(data.columns)]
    data = pd.concat([data, data_rev]).sort_index()

    data.reset_index().to_csv(output_csv, index=False)

@main.command()
@click.argument('input_csv')
@click.argument('output_csv')
def first_correct(input_csv, output_csv):
    data = pd.read_csv(input_csv)
    data = data.loc[data['rmsd1'] <= 2.05]
    data['native'] = data['rmsd2'] <= 2.05
    data['gscore'] = data['gscore2']
    data.reset_index().to_csv(output_csv, index=False)

@main.command()
@click.argument('input_csv')
@click.argument('output_csv')
def weight(input_csv, output_csv):
    data = pd.read_csv(input_csv)
    data = data.set_index(['protein', 'ligand1', 'ligand2'])
    
    # each ligand pair equal
    w = data.loc[data['native'] == 1, 'contact'].groupby(level=[0, 1, 2]).count()
    data.loc[data['native'] == 1, 'W_lig_pair'] = 1/w

    w = data.loc[data['native'] != 1, 'contact'].groupby(level=[0, 1, 2]).count()
    data.loc[data['native'] != 1, 'W_lig_pair'] = 1/w

    # each ligand equal
    w = data.reset_index().groupby(['protein', 'ligand2']).nunique()[['ligand1']]
    w = 1/w.rename(columns={'ligand1': 'W_ligand1'})
    data = data.join(w)

    # each protein equal (when multiplied by above)
    w = data.reset_index().groupby(['protein', 'ligand1']).nunique()[['ligand2']]
    w = 1/w.rename(columns={'ligand2': 'W_ligand2'})
    data = data.join(w)

    data.reset_index().to_csv(output_csv, index=False)

main()
