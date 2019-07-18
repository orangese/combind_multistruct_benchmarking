import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from shared_paths import proteins

with open('/scratch/PI/rondror/combind/mcss.csv', 'w') as fp:
	fp.write(','.join(['protein', 'ligand', 'rank', 'rmsd', 'gscore', 'mcss', 'mcss_atoms', 'l1_atoms', 'l2_atoms'])+'\n')
	for prot in proteins:
		print(prot)
		protein = Protein(prot)
		protein.lm.mcss.load_mcss()
		crystal_lig = "{}_crystal_lig".format(protein.lm.st)
		for ligand in protein.lm.get_xdocked_ligands(20):
			protein.docking[protein.lm.st].ligands = {}
			protein.docking[protein.lm.st].num_poses = {}
			try:
				protein.load_docking([crystal_lig], load_mcss=True, load_crystal = True)
				protein.load_docking([ligand], load_mcss=True, load_crystal = False)
			except:
				continue

			mcss_size = protein.lm.mcss.get_mcss_size(crystal_lig, ligand)

			mcss_atoms, l1_atoms, l2_atoms = protein.lm.mcss.get_mcss_and_ligand_sizes(crystal_lig, ligand)

			if mcss_size <= 0.5: continue

			for rank in range(min(100, protein.docking[protein.lm.st].num_poses[ligand])):
				pose = protein.docking[protein.lm.st].ligands[ligand].poses[rank]
				mcss = protein.lm.mcss.get_rmsd(crystal_lig, ligand, 0, rank)
				fp.write(','.join([str(x) for x in [prot, ligand, rank, pose.rmsd, pose.gscore, mcss, mcss_atoms, l1_atoms, l2_atoms]])+'\n')
