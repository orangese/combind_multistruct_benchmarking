from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pandas as pd

def get_fp(mol):
	return AllChem.GetMorganFingerprint(mol, 2)


import sys

ref = pd.read_csv('/oak/stanford/groups/rondror/users/jpaggi/DUDE/combind/{}/structures/pdb.csv'.format(sys.argv[1]))
assert len(ref) == 1
ref_smiles = ref.loc[0, 'SMILES']
ref_mol = Chem.MolFromSmiles(ref_smiles)
ref_fp = get_fp(ref_mol)

active = pd.read_csv('/oak/stanford/groups/rondror/users/jpaggi/DUDE/combind/{}/random_binders.csv'.format(sys.argv[1]))
active_fps = []
for i, ligand in active.iterrows():
	mol = Chem.MolFromSmiles(ligand['SMILES'])
	active_fps += [get_fp(mol)]

df = pd.read_csv('/oak/stanford/groups/rondror/users/jpaggi/DUDE/combind/{}/subset.csv'.format(sys.argv[1]))

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
		pass

print(df)
df.to_csv('/oak/stanford/groups/rondror/users/jpaggi/DUDE/combind/{}/subset_anno.csv'.format(sys.argv[1]), index=False)
