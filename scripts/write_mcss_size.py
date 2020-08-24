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
params = config.STATS['rd1_all']

proteins = utils.get_proteins(paths, [])

mcss_sizes = {}
for prot in proteins:
	print(prot)
	protein = Protein(prot, params, paths)
	protein.lm.mcss.load_mcss()
	crystal_lig = "{}_lig".format(protein.lm.st)
	for ligand in protein.lm.get_pdb(25):
		if ligand == crystal_lig: continue
		protein.lm.params['mcss_func'] = max
		mcss_sizes[ligand] = protein.lm.mcss.get_mcss_size(ligand, crystal_lig)

with open('stats_data/mcss_sizes.pkl', 'wb') as fp:
	pickle.dump(mcss_sizes, fp)
