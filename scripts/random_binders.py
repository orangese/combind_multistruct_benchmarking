import click
import pandas as pd
import numpy as np
import os
from subprocess import run
from glob import glob

@click.group()
def main():
	pass

@main.command()
@click.option('--n-train', default=10)
@click.option('--n-folds', default=10)
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
	print(input_csv)
	for cwd in glob(root + '/[0-9]'):
		if os.path.exists('{}/combind_preds.csv'.format(cwd)):
			continue
		run('$COMBINDHOME/main.py --data {} --ligands {}/binders.csv score {} all'.format(data, root, protein),
	    	shell=True, cwd=cwd)
		run('$COMBINDHOME/main.py --data {} --ligands {} screen {} all combind_preds.csv'.format(data, input_csv, protein),
	    	shell=True, cwd=cwd)

main()


"""
# Build and then test qsar models using Schrodinger's autoqsar tool

# Categorical model
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -build -i train.csv -y AFFINITY -cat -cuts 1000
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -test  -i test.csv  -pred autoqsar_preds.csv

# Quantitative model
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -build -i train.csv -y AFFINITY -num -log -scale 1.0e-9
$SCHRODINGER/utilities/autoqsar model.qzip -WAIT -test  -i test.csv  -pred autoqsar_preds.csv
"""