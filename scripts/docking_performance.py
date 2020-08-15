import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
import utils
import config

paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
         'DATA': '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020',
         'PDB': '{ROOT}/structures/pdb.csv'}
paths.update(config.PATHS)
paths = utils.resolve(paths)
params = config.STATS['shape']

proteins = utils.get_proteins(paths, [])

data1 = {}
for prot in proteins:
	protein = Protein(prot, params, paths)
	ligands = protein.lm.get_pdb()
	xtal = protein.lm.st + '_lig'
	protein.load_docking(ligands)
	for name, ligand in protein.docking.items():
		rmsds = [pose.rmsd for pose in ligand.poses]

		if rmsds:
			best = min(rmsds)
			top = rmsds[0]
		else:
			best = float('inf')
			top = float('inf')

		if name == xtal: continue

		data1[name] = (top, best)

		# if top <= 2.0:
		# 	print(prot, name, top, best)

params = config.STATS['rd1_all']
data2 = {}
for prot in proteins:
	protein = Protein(prot, params, paths)
	ligands = protein.lm.get_pdb()
	xtal = protein.lm.st + '_lig'
	protein.load_docking(ligands)
	for name, ligand in protein.docking.items():
		rmsds = [pose.rmsd for pose in ligand.poses]

		if rmsds:
			best = min(rmsds)
			top = rmsds[0]
		else:
			best = float('inf')
			top = float('inf')

		if name == xtal: continue

		data2[name] = (top, best)

		# if top <= 2.0:
		# 	print(prot, name, top, best)


for name in data1:
	print(name,
	      data1[name][0], data1[name][1],
	      data2[name][0], data2[name][1])
