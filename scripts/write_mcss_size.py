import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from settings import proteins, stats, paths

mcss_sizes = {}
for prot in proteins:
	print(prot)
	protein = Protein(prot, stats['stats41'], paths)
	protein.lm.mcss.load_mcss()
	crystal_lig = "{}_crystal_lig".format(protein.lm.st)
	for ligand in protein.lm.get_xdocked_ligands(20):
		mcss_sizes[ligand] = protein.lm.mcss.get_mcss_size(ligand, crystal_lig)

with open('/oak/stanford/groups/rondror/users/jpaggi/mcss_sizes.pkl', 'wb') as fp:
	pickle.dump(mcss_sizes, fp)
