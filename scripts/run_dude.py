"""
# Convert raw dud-e download into combind format.
python $COMBINDHOME/scripts/run_dude.py convert dude combind

# Prepare structure
combind structprep XTAL

# Dock screening libary
mkdir subset
mkdir subset/ligands
mkdir subset/docking
combind ligprep subset.smi subset/ligands
combind dock structures/grids/XTAL/XTAL.zip subset/docking subset/ligands/subset.maegz

for i in *; do cd $i; sbatch -p owners -t 24:00:00 --wrap="combind dock structures/grids/XTAL/XTAL.zip subset/docking subset/ligands/subset.maegz"; cd ..; done;

# Setup cross-validation sets.
mkdir $i/scores
mkdir $i/scores/rd1_all_15
python $COMBINDHOME/scripts/run_dude.py setup      $i/subset.smi $i/structures/pdb.smi $i/scores/rd1_all_15 --n-train 15 --n-folds 5

for i in *; do sbatch -p rondror --wrap="python $COMBINDHOME/scripts/run_dude.py shape $i/subset.smi $i/scores/rd1_all_15 $i" -J shape; done;
for i in *; do python $COMBINDHOME/scripts/run_dude.py similarity $i/subset.smi $i/scores/rd1_all_15 $i; done;

# Predict poses for selected binders

## Ligprep and docking
for i in */scores/rd1_all_3/*; do cd $i; mkdir bpp; cd -; done;
for i in */scores/rd1_all_15/[1-4]; do cd $i; sbatch -p owners --wrap="combind ligprep binder.smi bpp/ligands" -J ligprep; cd -; done;
for i in */scores/rd1_all_15/1; do cd $i; if [ ! -e bpp/binder.csv ]; then echo $i; sbatch  -p owners -t 06:00:00 --wrap="combind dock ../../../structures/grids/XTAL/XTAL.zip bpp/docking bpp/ligands/*/*.maegz" -J dock;  fi; cd -; done;
for i in *; do combind check-dock $i/scores/rd1_all_15/$j/binder.smi $i/structures/grids/XTAL/XTAL.zip $i/scores/rd1_all_15/$j/bpp; done;

for i in */scores/rd1_all_3/*; do cd $i; if [ ! -e bpp/docking/XTAL-to-XTAL/XTAL-to-XTAL_native_pv.maegz ]; then combind filter-native bpp/docking/XTAL-to-XTAL/XTAL-to-XTAL_pv.maegz ../../../structures/ligands/XTAL_lig.mae; fi; cd -; done;

## Featurization and scoring
for i in */scores/rd1_all_15/$j; do cd $i; combind featurize bpp bpp/docking/XTAL-to-XTAL/XTAL-to-XTAL_native_pv.maegz; cd -; done;

for i in */scores/rd1_all_15/[1-4]; do cd $i; if [ ! -e bpp/binder.csv ]; then sbatch  -p rondror -t 08:00:00 --wrap="combind featurize bpp bpp/docking/XTAL-to-XTAL/XTAL-to-XTAL_native_pv.maegz bpp/docking/CHEMBL*/*pv.maegz" -J features; fi; cd -; done;
for i in */scores/rd1_all_15/*; do cd $i; if [ ! -e bpp/binder.csv ]; then echo $i; sbatch -p rondror -J score --wrap="combind pose-prediction bpp bpp/binder.csv --xtal XTAL-to-XTAL_native_pv --gc50 -8.0 --alpha 1.0 --features shape,mcss,hbond,saltbridge,contact"; fi; cd -; done;
for i in */scores/rd1_all_15/*; do cd $i; if [ ! -e bpp/binder_pv.maegz ]; then echo $i; combind extract-top-poses bpp/binder.csv bpp/docking; fi; cd -; done;

# ComBind screening
for i in */scores/rd1_all_15/*; do cd $i; mkdir screen; cd -; done;
for i in */scores/rd1_all_15/1; do cd $i; echo $i; if [ ! -e combind.csv ]; then sbatch --begin=now+0hours -p rondror -t 12:00:00 --mem=32GB --wrap="combind featurize screen ../../../subset/docking/subset-to-XTAL/subset-to-XTAL_pv.maegz bpp/binder_pv.maegz --no-mcss --screen" -J screen; fi; cd -; done;

for i in */scores/rd1_all_15/*; do cd $i; if [ ! -e combind.csv ]; then echo $i; python ~/combind/scripts/run_dude.py combind-screen; fi; cd -; done;

for i in */scores/rd1_all_15/*; do cd $i; if [ ! -e combind.csv ]; then echo $i; sbatch -p rondror --wrap="python ~/combind/scripts/run_dude.py combind-screen" -J apply; fi; cd -; done;

for i in */scores/rd1_all_*/0; do cd $i; if [ ! -e combind_dyn.csv ]; then echo $i; sbatch -p rondror --wrap="python ~/combind/scripts/run_dude.py combind-screen-dyn" -J apply; fi; cd -; done;

for k in 0 1 3 5 10; do for j in {0..4}; do echo $k $j $(ls */scores/*_$k/$j/sim*.csv | wc -l); done; done;

util.cbay; util.cbac *.CHEMBL*; util.cbam binder*.*;
as cartoon; show lines; show sticks, het; hide everything, element H and not (element O+N+S extend 1)

reni/scores/rd1_all_15/2
reni/scores/rd1_all_15/3
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

@click.group()
def main():
    pass

@main.command()
@click.argument('dude_dir', type=click.Path(exists=True))
@click.argument('combind_dir', type=click.Path(exists=True))
@click.argument('xtal_csv', type=click.Path(exists=True))
def convert(dude_dir, combind_dir, xtal_csv):
    xtal_df = pd.read_csv(xtal_csv)
    xtal_df = xtal_df.set_index('Entry')
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

        actives.loc[actives.CHEMBL.isna(), 'CHEMBL'] = \
            ['CHEMBLX'+str(x) for x in actives.loc[actives.CHEMBL.isna(), 'ID']]
        actives['ID'] = actives['CHEMBL']


        actives = actives.loc[:, ['SMILES', 'ID']]
        decoys = decoys.loc[:, ['SMILES', 'ID']]
        
        all_ligands = pd.concat([actives, decoys])

        np.random.seed(42)
        subset_ligands = pd.concat([actives, decoys.sample(1000)])

        with StructureReader(dude_path + '/receptor.pdb') as st:
            receptor = list(st)[0]

        with StructureReader(dude_path + '/crystal_ligand.mol2') as st:
            ligand = list(st)[0]

        os.mkdir(combind_path)
        os.mkdir(combind_path + '/structures')
        os.mkdir(combind_path + '/structures/raw')

        all_ligands.to_csv(combind_path + '/all.smi', index=False, sep=' ')
        subset_ligands.to_csv(combind_path + '/subset.smi', index=False, sep=' ')

        ligand.write(combind_path + '/structures/raw/XTAL_lig.mae')
        receptor.write(combind_path + '/structures/raw/XTAL_prot.mae')

        with open(combind_path + '/structures/pdb.smi', 'w') as fp:
            fp.write('ID SMILES\n')
            fp.write('XTAL {}\n'.format(xtal_df.loc[protein.upper(), 'SMILES']))

@main.command()
@click.option('--n-train', default=10)
@click.option('--n-folds', default=7)
@click.argument('input_csv')
@click.argument('xtal')
@click.argument('root')
def setup(input_csv, xtal, root, n_train, n_folds):
    np.random.seed(42)
    df = pd.read_csv(input_csv, sep=' ')
    xtal = pd.read_csv(xtal, sep=' ').loc[:1]
    
    for i in range(n_folds):
        cwd = '{}/{}'.format(root, i)
        binder_smi = '{}/{}/binder.smi'.format(root, i)
        
        if os.path.exists(cwd):
            print(cwd, 'exists. not overwriting.')
            continue

        binders = df.loc[df['ID'].str.contains('CHEMBL')]
        binders = binders.sample(n_train)
        binders = pd.concat([binders, xtal])

        os.mkdir(cwd)
        binders.to_csv(binder_smi, index=False, sep=' ')

@main.command()
def combind_screen():
    if not os.path.exists('screen/screen.npy'):
        print('score')
        run('combind screen screen/screen.npy screen/gscore/subset-to-XTAL_pv.npy --ifp-fname screen/ifp-pair/{}-subset-to-XTAL_pv-and-binder_pv.npy --shape-fname screen/shape/shape-subset-to-XTAL_pv-and-binder_pv.npy', shell=True)

    if os.path.exists('screen/screen.npy') and not os.path.exists('screen/screen_pv.maegz'):
        print('apply')
        run('combind apply-scores ../../../subset/docking/subset-to-XTAL/subset-to-XTAL_pv.maegz screen/screen.npy screen/screen_pv.maegz', shell=True)

    if os.path.exists('screen/screen_pv.maegz') and not os.path.exists('screen/screen_combind_pv.maegz'):
        print('extract-combind')
        run('$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score  -o screen/screen_combind_pv.maegz screen/screen_pv.maegz', shell=True)

    if os.path.exists('screen/screen_pv.maegz') and not os.path.exists('screen/screen_glide_pv.maegz'):
        print('extract-glide')
        run('$SCHRODINGER/utilities/glide_sort -best_by_title -o screen/screen_glide_pv.maegz screen/screen_pv.maegz', shell=True)

    if os.path.exists('screen/screen_combind_pv.maegz') and not os.path.exists('combind.csv'):
        print('csv-combind')
        run('combind scores-to-csv screen/screen_combind_pv.maegz combind.csv', shell=True)

    if os.path.exists('screen/screen_glide_pv.maegz') and not os.path.exists('glide.csv'):
        print('csv-glide')
        run('combind scores-to-csv screen/screen_glide_pv.maegz glide.csv', shell=True)

@main.command()
def combind_screen_dyn():
    if not os.path.exists('screen/screen_dyn.npy'):
        print('score')
        run('combind screen screen/screen_dyn.npy screen/gscore/subset-to-XTAL_pv.npy --ifp-fname screen/ifp-pair/{}-subset-to-XTAL_pv-and-binder_pv.npy --shape-fname screen/shape/shape-subset-to-XTAL_pv-and-binder_pv.npy', shell=True)

    if os.path.exists('screen/screen_dyn.npy') and not os.path.exists('screen/screen_dyn_pv.maegz'):
        print('apply')
        run('combind apply-scores ../../../subset/docking/subset-to-XTAL/subset-to-XTAL_pv.maegz screen/screen_dyn.npy screen/screen_dyn_pv.maegz', shell=True)

    if os.path.exists('screen/screen_dyn_pv.maegz') and not os.path.exists('screen/screen_dyn_combind_pv.maegz'):
        print('extract-combind')
        run('$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score  -o screen/screen_dyn_combind_pv.maegz screen/screen_dyn_pv.maegz', shell=True)

    if os.path.exists('screen/screen_dyn_combind_pv.maegz') and not os.path.exists('combind_dyn.csv'):
        print('csv-combind')
        run('combind scores-to-csv screen/screen_dyn_combind_pv.maegz combind_dyn.csv', shell=True)


@main.command()
def combind_screen_prob():
    if not os.path.exists('screen/screen_prob.npy'):
        print('score')
        run('combind screen screen/screen_prob.npy screen/gscore/subset-to-XTAL_pv.npy --ifp-fname screen/ifp-pair/{}-subset-to-XTAL_pv-and-binder_pv.npy --shape-fname screen/shape/shape-subset-to-XTAL_pv-and-binder_pv.npy --pose-prob-fname bpp/binder.csv', shell=True)

    if os.path.exists('screen/screen_prob.npy') and not os.path.exists('screen/screen_prob_pv.maegz'):
        print('apply')
        run('combind apply-scores ../../../subset/docking/subset-to-XTAL/subset-to-XTAL_pv.maegz screen/screen_prob.npy screen/screen_prob_pv.maegz', shell=True)

    if os.path.exists('screen/screen_prob_pv.maegz') and not os.path.exists('screen/screen_prob_combind_pv.maegz'):
        print('extract-combind')
        run('$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score  -o screen/screen_prob_combind_pv.maegz screen/screen_prob_pv.maegz', shell=True)

    if os.path.exists('screen/screen_prob_combind_pv.maegz') and not os.path.exists('combind_prob.csv'):
        print('csv-combind')
        run('combind scores-to-csv screen/screen_prob_combind_pv.maegz combind_prob.csv', shell=True)

def get_fp(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

@main.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
@click.argument('data_root')
def similarity(input_csv, root, data_root):
    ref = pd.read_csv('{}/structures/pdb.smi'.format(data_root), sep=' ')
    assert len(ref) == 1
    ref_smiles = ref.loc[0, 'SMILES']
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_fp = get_fp(ref_mol)

    df = pd.read_csv(input_csv, sep=' ')
    for cwd in glob(root + '/[0-9]'):
        sim_fname = '{}/similarity.csv'.format(cwd)
        if os.path.exists(sim_fname):
            print(sim_fname, 'exists. not overwriting.')
            continue

        active = pd.read_csv('{}/binder.smi'.format(cwd), sep=' ')
        active_fps = []
        for i, ligand in active.iterrows():
            mol = Chem.MolFromSmiles(ligand['SMILES'])
            active_fps += [get_fp(mol)]

        df['XTAL_sim'] = 0
        df['active_sim_max'] = 0
        df['active_sim_mean'] = 0
        for i, ligand in df.iterrows():
            try:
                mol = Chem.MolFromSmiles(ligand['SMILES'])
                fp = get_fp(mol)
                df.loc[i, 'XTAL_sim'] = DataStructs.TanimotoSimilarity(fp, ref_fp)
                active_sims = [DataStructs.TanimotoSimilarity(fp, active_fp)
                               for active_fp in active_fps]
                df.loc[i, 'active_sim_max'] = max(active_sims)
                df.loc[i, 'active_sim_mean'] = np.mean(active_sims)
            except:
                pass
        df.to_csv(sim_fname, index=False)

def extract(root, ligands, out):
    paths = ['{}/ligands/{}/{}.maegz'.format(root, ligand, ligand)
             for ligand in ligands]

    cmd = 'python $COMBINDHOME/scripts/shape_screen.py extract --best {} {}'.format(out, ' '.join(paths))
    os.system(cmd)

@main.command()
@click.argument('input_csv', type=click.Path(exists=True))
@click.argument('root')
@click.argument('data_root')
def shape(input_csv, root, data_root):
    for cwd in glob(root + '/[0-9]'):
        template = os.path.abspath('{}/structures/ligands/XTAL_lig.mae'.format(data_root))
        test_mae = os.path.abspath('{}/subset/ligands/subset.maegz'.format(data_root))

        active_smi = os.path.abspath('{}/binder.smi'.format(cwd))
        active_mae = os.path.abspath('{}/shape/binder.maegz'.format(cwd))
        active_align_mae = os.path.abspath('{}/shape/binder-to-XTAL_lig_align.maegz'.format(cwd))
        test_align_csv = os.path.abspath('{}/shape/subset-to-binder-to-XTAL_lig_align_align.csv'.format(cwd))
        shape_csv = '{}/shape.csv'.format(cwd)

        if os.path.exists(shape_csv):
            continue

        if os.path.exists('{}/shape'.format(cwd)):
            print('Shape not complete but temp directory exists. '
                  'Remove directory and try again.')
            print('{}/shape'.format(cwd))
            continue

        os.mkdir('{}/shape'.format(cwd))

        active = pd.read_csv(active_smi, sep=' ')['ID']
        extract('{}/bpp'.format(cwd), active, active_mae)

        cmd = 'python $COMBINDHOME/scripts/shape_screen.py screen {} {}'.format(template, active_mae)
        run(cmd, cwd=cwd+'/shape', shell=True)

        cmd = 'python $COMBINDHOME/scripts/shape_screen.py screen {} {}'.format(active_align_mae, test_mae)
        run(cmd, cwd=cwd+'/shape', shell=True)

        df = pd.read_csv(test_align_csv)
        df = df.set_index('ID')
        mean = df.groupby('ID')['score'].mean()
        df = df.groupby('ID')[['score']].max()
        df['SHAPE_max'] = df.score
        df['SHAPE_mean'] = mean
        df = df[['SHAPE_mean', 'SHAPE_max']]
        df.to_csv(shape_csv)

main()
