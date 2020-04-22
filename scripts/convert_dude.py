import pandas as pd
import numpy as np
import click
from glob import glob
import os
from schrodinger.structure import StructureReader

@click.command()
@click.argument('dude_dir', type=click.Path(exists=True))
@click.argument('combind_dir', type=click.Path(exists=True))
def main(dude_dir, combind_dir):
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

main()