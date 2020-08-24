import numpy as np
from containers import Protein
import config
import utils
from score.pairs import LigPair
from glob import glob
import os
import pandas as pd


def load_ligand_pair(prot, protein, ligand1, ligand2):
    lig_pair = LigPair(prot.docking[ligand1],
                       prot.docking[ligand2],
                       _interactions,
                       prot.lm.mcss  if 'mcss' in interactions else None,
                       prot.lm.shape if 'shape' in interactions else None,
                       settings['max_poses'])

    n1 = min(len(lig_pair.l1.poses), settings['max_poses'])
    n2 = min(len(lig_pair.l2.poses), settings['max_poses'])

    X = []
    for r1 in range(n1):
        for r2 in range(n2):
        	features = [lig_pair.get_feature(interaction, r1, r2)
        		        for interaction in interactions]
        	X += [[prot.protein,
        		   ligand1, ligand2,
        		   prot.docking[ligand1].poses[r1].gscore,
        		   prot.docking[ligand2].poses[r2].gscore,
        		   prot.docking[ligand1].poses[r1].rmsd,
        		   prot.docking[ligand2].poses[r2].rmsd,
        		   r1, r2]+features]
    X = pd.DataFrame(X, columns=['protein',
                     'ligand1', 'ligand2',
                     'gscore1', 'gscore2',
                     'rmsd1', 'rmsd2',
                     'rank1', 'rank2']+interactions)
    return X

def load_protein(protein):
	prot = Protein(protein, settings, paths)
	ligands = prot.lm.get_xdocked_ligands(10)
	print(ligands)
	prot.load_docking(ligands, load_fp=True,
	                  load_mcss='mcss' in interactions,
	                  load_shape='shape' in interactions)
	df = []
	for j, ligand1 in enumerate(ligands):
		for ligand2 in ligands[j+1:]:
			df += [load_ligand_pair(prot, protein, ligand1, ligand2)]
	return pd.concat(df)

stats_version = 'rd1'
data = '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems'
ligands = '{ROOT}/structures/pdb.csv'

settings = config.STATS[stats_version]
paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
		 'DATA': data,
		 'PDB': ligands}
paths.update(config.PATHS)
paths = utils.resolve(paths)

_interactions = config.FEATURE_DEFS
interactions = list(_interactions)

for protein in utils.get_proteins(paths, []):
	if os.path.exists('pairs/{}.csv'.format(protein)): continue
	print(protein)
	load_protein(protein).to_csv('pairs/{}.csv'.format(protein))
