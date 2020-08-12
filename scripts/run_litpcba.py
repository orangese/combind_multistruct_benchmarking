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
@click.argument('raw_dir', type=click.Path(exists=True))
@click.argument('combind_dir', type=click.Path(exists=True))
def convert(raw_dir, combind_dir):
    for raw_path in glob(raw_dir + '/*'):
        protein = raw_path.split('/')[-1]
        combind_path = combind_dir + '/' + protein

        if os.path.exists(combind_path):
            print(combind_path, 'already exists')
            continue

        actives = pd.read_csv(raw_path + '/actives.smi',
                              sep=' ', names=['SMILES', 'ID'])
        decoys = pd.read_csv(raw_path + '/inactives.smi',
                              sep=' ', names=['SMILES', 'ID'])

        actives['AFFINITY'] = 1
        decoys['AFFINITY'] = 1e6

        all_ligands = pd.concat([actives, decoys])

        np.random.seed(42)
        subset_ligands = pd.concat([actives.sample(1000), decoys.sample(1000)])

        with StructureReader(sorted(glob(raw_path + '/*_ligand.mol2'))[0]) as st:
            ligand = list(st)[0]

        with StructureReader(sorted(glob(raw_path + '/*_protein.mol2'))[0]) as st:
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

main()
