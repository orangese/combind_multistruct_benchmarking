"""
python scripts/pick_helpers.py /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/structures/pdb.csv /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/chembl/CHEMBL*.csv /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/chembl --criteria affinity
"""

import subprocess
from glob import glob
import pandas as pd
import os

root = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020'
path = '{root}/{protein}/chembl'
chembl = 'stats_data/systems.txt'

cmd = 'sbatch -p rondror -t 02:00:00 --wrap="python ~/combind/scripts/chembl.py {chembl}"'

systems = pd.read_csv(chembl, sep='\t')
for _, system in systems.iterrows():
	protein = system['UNIPROT']
	chembl_id = system['CHEMBL']

	cwd = path.format(root=root, protein=protein)

	if not os.path.exists(cwd):
		os.mkdir(cwd)
	
	if glob(cwd + '/CHEMBL*'):
		continue
	
	if protein == 'Q05586-Q12879':
		_cmd = cmd.format(chembl=chembl_id + ' --protein-complex')
	else:
		_cmd = cmd.format(chembl=chembl_id)
	print(cwd, _cmd)
	subprocess.run(_cmd, shell=True, cwd=cwd)


systems = pd.read_csv(chembl, sep='\t')
for _, system in systems.iterrows():
	protein = system['UNIPROT']
	chembl_id = system['CHEMBL']
	cwd = path.format(root=root, protein=protein)

	csv = glob(cwd + '/CHEMBL*_all.csv')
	if not csv: continue
	assert len(csv) == 1, csv
	csv = csv[0]

	df = pd.read_csv(csv)
	df['AFFINITY'] = df['standard_value']
	df['SMILES'] = df['canonical_smiles']
	df['ID'] = df['ligand_chembl_id']

	a = df['standard_value']
	df[a < 10**3].to_csv(csv.replace('.csv', '_nM.csv'), index=False)
	df[(a >= 10**3) & (a < 10**4)].to_csv(csv.replace('.csv', '_uM.csv'))
	df[a >= 10**4].to_csv(csv.replace('.csv', '_mM.csv'), index=False)


cmd = 'python scripts/pick_helpers.py {root}/structures/pdb.csv {root}/chembl/CHEMBL*_nM.csv {root}/chembl --criteria affinity'
systems = pd.read_csv(chembl, sep='\t')
for _, system in systems.iterrows():
	protein = system['UNIPROT']
	chembl_id = system['CHEMBL']
	_root = '{}/{}'.format(root, protein)
	_cmd = cmd.format(root=_root)
	print(_cmd)
	subprocess.run(_cmd, shell=True)


