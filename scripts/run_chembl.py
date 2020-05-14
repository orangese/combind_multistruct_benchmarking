"""
# Start with structural data in combind format. Here I used the proteins
# I had already studied for the pose prediction paper.
# Additionally, you should write the chembl target id that you wish to
# use for each target to a file. It's hard to do this in an automated way
# as sometimes you want to use chembl data from a different organism than
# the structural data.

# Query CHEMBL database
python $COMBINDHOME/scripts/run_chembl.py get-chembl

# Split the chembl data into IC50 and Ki. Make a link to whichever one
# has more data under "chosen.csv"
# This also fixes some column names and creates a LOGAFFINITY column
python $COMBINDHOME/scripts/run_chembl.py process-chembl $i/chembl

# Dock ligands and compute interaction fingerprints
./main.py --ligands '{ROOT}/chembl/chosen.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL prepare prep-structs
./main.py --ligands '{ROOT}/chembl/chosen.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL prepare prep-ligands
./main.py --ligands '{ROOT}/chembl/chosen.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL prepare dock
./main.py --ligands '{ROOT}/chembl/chosen.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL prepare ifp

# Setup cross-validation sets.
mkdir $i/scores
mkdir $i/scores/rd1
python $COMBINDHOME/scripts/run_chembl.py setup $i/chembl/chosen.csv $i/scores/rd1

# Compute MCSS for each set of "helper ligands".
./main.py --ligands '{ROOT}/scores/rd1/0/binder.csv' --data /oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL prepare mcss

# Make predictions with ComBind
python $COMBINDHOME/scripts/run_chembl.py combind $i/chembl/chosen.csv $i/scores/rd1 . $i

# Make predictions with AUTOQSAR
python $COMBINDHOME/scripts/run_chembl.py autoqsar $i/chembl/chosen.csv $i/scores/rd1

# Annotate with similarity to the XTAL ligand and helper ligands.
python $COMBINDHOME/scripts/run_chembl.py similarity $i/chembl/chosen.csv $i/scores/rd1 . $i

# Useful for submitting jobs for each protein
sbatch -p rondror -t 03:00:00 -J $i -o ~/temp/$i.log --wrap="$CMD"
for i in $(ls --color=none); do $CMD; done;
"""

import os
from glob import glob
import click
import numpy as np
import pandas as pd
from subprocess import run
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

from rdkit.Chem.rdmolfiles import MaeMolSupplier

@click.group()
def main():
    pass

@main.command()
def get_chembl():
    with open('scripts/chembl/chembl.txt') as fp:
        for line in fp:
            protein, chembl = line.strip().split()

            wd = '/oak/stanford/groups/rondror/users/jpaggi/VS/CHEMBL/{}/chembl'.format(protein)
            cmd = 'python /home/users/jpaggi/combind_dev/scripts/chembl/chembl.py {}'.format(wd, chembl)
            if chembl == 'CHEMBL1907604':
                cmd += ' --protein_complex'
            else:
                continue

            cmd = 'sbatch -p rondror -t 01:00:00 -J chembl -D {} --wrap="{}"'.format(wd, cmd)

            if not os.path.exists(wd):
                os.mkdir(wd)

            exists = glob(wd + '/CHEMBL*')

            if not exists:
                os.system(cmd)

@main.command()
@click.argument('chembl_root')
def process_chembl(chembl_root):
    chembl_fname = glob(chembl_root + '/CHEMBL*.csv')
    if not len(chembl_fname):
        print('No CHEMBL file.')
        return
    if len(chembl_fname) > 1:
        print('Multiple chembl files')
        return
    chembl_fname = chembl_fname[0]

    df = pd.read_csv(chembl_fname)
    df = df.rename(columns={'ligand_chembl_id': 'ID',
                            'standard_value': 'AFFINITY',
                            'canonical_smiles': 'SMILES'})

    df['LOGAFFINITY'] = np.log10(df['AFFINITY']) - 9

    cols = df.columns.tolist()
    cols.remove('ID')
    cols.remove('LOGAFFINITY')
    cols.remove('AFFINITY')
    cols.remove('SMILES')
    cols = ['SMILES', 'ID', 'AFFINITY', 'LOGAFFINITY'] + cols
    df = df[cols]

    Ki = df.loc[df['standard_type'] == 'Ki']
    IC50 = df.loc[df['standard_type'] == 'IC50']

    Ki.to_csv(chembl_root + '/Ki.csv', index=False)
    IC50.to_csv(chembl_root + '/IC50.csv', index=False)

    if Ki.shape[0] > IC50.shape[0]:
        os.system('ln -s Ki.csv {}/chosen.csv'.format(chembl_root))
    else:
         os.system('ln -s IC50.csv {}/chosen.csv'.format(chembl_root))


@main.command()
@click.option('--n-train', default=10)
@click.option('--n-folds', default=7)
@click.argument('input_csv')
@click.argument('root')
def setup(input_csv, root, n_train, n_folds):
    np.random.seed(42)
    df = pd.read_csv(input_csv)
    
    for i in [-7, -6]:
        if sum(df['LOGAFFINITY'] <= i) > 30:
            affinity_cut = i
            print(affinity_cut)
            break
    else:
        print('exit.')
        return

    for i in range(n_folds):
        cwd = '{}/{}'.format(root, i)
        if os.path.exists(cwd):
            print(cwd, 'exists. not overwriting.')
            continue
        
        binder_csv = '{}/{}/binder.csv'.format(root, i)
        binders = df.loc[df['LOGAFFINITY'] < affinity_cut]
        binders = binders.sample(n_train)
        os.mkdir(cwd)
        binders.to_csv(binder_csv, index=False)

@main.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
def autoqsar(input_csv, root):
    input_csv = os.path.abspath(input_csv)
    root = os.path.abspath(root)
    for cwd in glob(root + '/[0-9]'):
        if os.path.exists('{}/autoqsar_preds.csv'.format(cwd)):
            continue
        run('$SCHRODINGER/utilities/autoqsar autoqsar.qzip -WAIT -build -i binder.csv -y AFFINITY -num -log -scale 1.0e-9',
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
        if not os.path.exists('{}/combind_preds.csv'.format(cwd)):
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
    ref_pdb = glob('{}/{}/docking/grids/*'.format(data, protein))
    assert len(ref_pdb) == 1
    ref_pdb = ref_pdb[0].split('/')[-1]
    ref_fname = '{}/{}/structures/ligands/{}_lig.mae'.format(data, protein, ref_pdb)
    print(ref_pdb, ref_fname)
    ref_mol = list(MaeMolSupplier(ref_fname))
    assert len(ref_mol) == 1
    ref_fp = get_fp(ref_mol[0])
    
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
