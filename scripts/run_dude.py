"""
# Convert raw dud-e download into combind format.
cd /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE
python $COMBINDHOME/scripts/run_dude.py convert dude/dud38 combind
python $COMBINDHOME/scripts/run_dude.py convert dude/diverse combind
python $COMBINDHOME/scripts/run_dude.py convert dude/gpcr combind

# Dock ligands and compute interaction fingerprints
./main.py --ligands '{ROOT}/subset.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE/combind prepare prep-structs
./main.py --ligands '{ROOT}/subset.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE/combind prepare prep-ligands
./main.py --ligands '{ROOT}/subset.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE/combind prepare dock
./main.py --ligands '{ROOT}/subset.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE/combind prepare ifp

# Setup cross-validation sets.
mkdir $i/scores
mkdir $i/scores/subset10_rd1
python $COMBINDHOME/scripts/run_dude.py setup $i/subset.csv $i/scores/subset10_rd1

# Compute MCSS for each set of "helper ligands".
./main.py --ligands '{ROOT}/scores/subset10_rd1/0/binder.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/DUDE/combind prepare mcss

# Make predictions with ComBind
python $COMBINDHOME/scripts/run_dude.py combind $i/subset.csv $i/scores/subset10_rd1 . $i

# Make predictions with AUTOQSAR
python $COMBINDHOME/scripts/run_dude.py autoqsar $i/subset.csv $i/scores/subset10_rd1

# Annotate with similarity to the XTAL ligand and helper ligands.
python $COMBINDHOME/scripts/run_dude.py similarity $i/subset.csv $i/scores/subset10_rd1 . $i

# Useful for submitting jobs for each protein
sbatch -p rondror -t 03:00:00 -J $i -o ~/temp/$i.log --wrap="$CMD"
for i in $(ls --color=none); do $CMD; done;
"""

import os
import sys
import click
import pandas as pd
import numpy as np
from glob import glob
from subprocess import run

import scipy.stats
import sklearn.metrics

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import generate_smiles

@click.group()
def main():
    pass

@main.command()
@click.argument('dude_dir', type=click.Path(exists=True))
@click.argument('combind_dir', type=click.Path(exists=True))
def convert(dude_dir, combind_dir):
    for dude_path in glob(dude_dir + '/*'):
        protein = dude_path.split('/')[-1]
        combind_path = combind_dir + '/' + protein

        if os.path.exists(combind_path):
            print(combind_path, 'already exists')
            continue

        actives = pd.read_csv(dude_path + '/' + 'actives_final.ism',
                              sep=' ', names=['SMILES', 'ID', 'CHEMBL'])
        decoys = pd.read_csv(dude_path + '/' + 'decoys_final.ism',
                              sep=' ', names=['SMILES', 'ID', 'CHEMBL'])

        actives['AFFINITY'] = 1
        decoys['AFFINITY'] = 1e6

        all_ligands = pd.concat([actives, decoys])

        np.random.seed(42)
        subset_ligands = pd.concat([actives, decoys.sample(1000)])

        with StructureReader(dude_path + '/crystal_ligand.mol2') as st:
            ligand = list(st)[0]

        with StructureReader(dude_path + '/receptor.pdb') as st:
            receptor = list(st)[0]

        os.mkdir(combind_path)
        os.mkdir(combind_path + '/structures')
        os.mkdir(combind_path + '/structures/raw')

        all_ligands.to_csv(combind_path + '/all.csv', index=False)
        subset_ligands.to_csv(combind_path + '/subset.csv', index=False)

        ligand.write(combind_path + '/structures/raw/XTAL_lig.mae')
        receptor.write(combind_path + '/structures/raw/XTAL_prot.mae')

        with open(combind_path + '/structures/pdb.csv', 'w') as fp:
            fp.write('ID,SMILES\n')
            fp.write('XTAL,{}\n'.format(generate_smiles(ligand)))

@main.command()
@click.option('--n-train', default=10)
@click.option('--n-folds', default=7)
@click.option('--affinity-cut', default=1000)
@click.argument('input_csv')
@click.argument('root')
def setup(input_csv, root, n_train, n_folds, affinity_cut):
    np.random.seed(42)
    df = pd.read_csv(input_csv)
    for i in range(n_folds):
        cwd = '{}/{}'.format(root, i)
        if os.path.exists(cwd):
            print(cwd, 'exists. not overwriting.')
            continue

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
        if not os.path.exists('{}/combind_poses.sc'.format(cwd)):
            run('$COMBINDHOME/main.py --data {} --ligands {}/binder.csv score {} all --pose-fname combind_poses.sc'.format(data, cwd, protein),
                shell=True, cwd=cwd)

        if not os.path.exists('{}/combind_preds.csv'.format(cwd)):
            run('$COMBINDHOME/main.py --data {} --ligands {} screen {} all --pose-fname combind_poses.sc --score-fname combind_preds.csv'.format(data, input_csv, protein),
                shell=True, cwd=cwd)

def get_fp(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

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
        sim_fname = '{}/similarity.csv'.format(cwd)
        if os.path.exists(sim_fname):
            print(sim_fname, 'exists. not overwriting.')
            continue

        active = pd.read_csv('{}/binder.csv'.format(cwd))
        active_fps = []
        for i, ligand in active.iterrows():
            mol = Chem.MolFromSmiles(ligand['SMILES'])
            active_fps += [get_fp(mol)]

        decoy = pd.read_csv('{}/decoy.csv'.format(cwd))
        decoy_fps = []
        for i, ligand in decoy.iterrows():
            mol = Chem.MolFromSmiles(ligand['SMILES'])
            decoy_fps += [get_fp(mol)]

        df['XTAL_sim'] = 0
        df['active_sim'] = 0
        df['decoy_sim'] = 0
        for i, ligand in df.iterrows():
            try:
                mol = Chem.MolFromSmiles(ligand['SMILES'])
                fp = get_fp(mol)
                df.loc[i, 'XTAL_sim'] = DataStructs.TanimotoSimilarity(fp, ref_fp)
                df.loc[i, 'active_sim'] = max(DataStructs.TanimotoSimilarity(fp, active_fp)
                                              for active_fp in active_fps)
                df.loc[i, 'decoy_sim'] = max(DataStructs.TanimotoSimilarity(fp, decoy_fp)
                                             for decoy_fp in decoy_fps)
            except:
                pass
        df.to_csv(sim_fname, index=False)
main()
