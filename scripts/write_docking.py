import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from settings import stats, proteins, paths

fname = '/oak/stanford/groups/rondror/users/jpaggi/docking.csv'

with open(fname, 'w') as fp:
	fp.write(','.join(['protein', 'ligand', 'rank', 'rmsd','gscore', 'emodel'])+'\n')
	for prot in proteins:
		print(prot)
		protein = Protein(prot, stats['stats41'], paths)
		ligands = protein.lm.get_xdocked_ligands(20)
		protein.load_docking(ligands)
		docking = protein.docking[protein.lm.st]

		for lig in ligands:
			for rank, pose in enumerate(docking[lig].poses):
				pose_info = [prot, lig, rank,
						     pose.rmsd, pose.gscore, pose.emodel]
				fp.write(','.join([str(x) for x in pose_info])+'\n')
