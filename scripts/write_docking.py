import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from shared_paths import shared_paths, proteins

fname = '/scratch/PI/rondror/combind/docking/{}.csv'.format(shared_paths['docking'])

with open(fname, 'w') as fp:
	fp.write(','.join(['protein', 'ligand', 'rank', 'rmsd','gscore', 'emodel'])+'\n')
	for prot in proteins:
		print(prot)
		protein = Protein(prot)
		protein.lm.mcss.load_mcss()
		ligands = protein.lm.get_xdocked_ligands(20)
		protein.load_docking(ligands)
		docking = protein.docking[protein.lm.st]

		for lig in ligands:
			for rank, pose in enumerate(docking.ligands[lig].poses):
				pose_info = [prot, lig, rank,
						     pose.rmsd, pose.gscore, pose.emodel]
				fp.write(','.join([str(x) for x in pose_info])+'\n')
