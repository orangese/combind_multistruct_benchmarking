import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from score.pairs import LigPair

prot = sys.argv[1]
fname = '/scratch/PI/rondror/combind/features/{}.csv'.format(prot, 'w')

features = ['mcss', 'hbond', 'sb', 'contact', 'pipi']
max_poses = 20

with open(fname, 'w') as fp:
	fp.write(','.join(['protein', 'ligand1', 'ligand2', 'rank1', 'rank2',
	                  'rmsd1', 'rmsd2', 'gscore1', 'gscore2']
	                  + features)+'\n')
	protein = Protein(prot)
	protein.lm.mcss.load_mcss()
	ligands = protein.lm.get_xdocked_ligands(20)
	protein.load_docking(ligands, load_mcss=True, load_fp=True)

	for i, lig1 in enumerate(ligands):
		for lig2 in ligands[i+1:]:
			lp = LigPair(protein.docking[protein.lm.st].ligands[lig1],
			             protein.docking[protein.lm.st].ligands[lig2],
			             features, protein.lm.mcss, max_poses)
			for (r1, r2), pp in lp.pose_pairs.items():
				pp_features = [lp.get_feature(f, r1, r2) for f in features]
				pose_info = [prot,
							 lig1, lig2,
							 r1, r2,
							 pp.pose1.rmsd, pp.pose2.rmsd,
							 pp.pose1.gscore, pp.pose2.gscore]
				fp.write(','.join([str(x) for x in pose_info+pp_features])+'\n')
