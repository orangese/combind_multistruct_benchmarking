import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from settings import proteins, stats, paths

version = sys.argv[1]

mcss_sizes = {}
for prot in proteins:
	print(prot)
	protein = Protein(prot, stats[version], paths)
	protein.lm.mcss.load_mcss()
	crystal_lig = "{}_crystal_lig".format(protein.lm.st)
	for ligand in protein.lm.get_xdocked_ligands(20):
		mcss_sizes[ligand] = protein.lm.mcss.get_mcss_size(ligand, crystal_lig)

with open('/oak/stanford/groups/rondror/users/jpaggi/{}_mcss_sizes.pkl'.format(version), 'wb') as fp:
	pickle.dump(mcss_sizes, fp)
