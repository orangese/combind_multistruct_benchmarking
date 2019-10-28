from glob import glob
import sys
import os

sys.path.append(os.getenv('COMBINDHOME'))
from settings import stats, proteins, paths
from containers import Protein


total, mcss = 0, 0
for protein in proteins:
	print(protein)
	prot = Protein(protein, stats['stats41'], paths)
	ligands = prot.lm.get_xdocked_ligands(20)
	prot.lm.mcss.load_mcss()

	for i, ligand1 in enumerate(ligands):
		for ligand2 in ligands[i+1:]:
			total += 1
			m = prot.lm.mcss.MCSSs['{}-{}'.format(ligand1, ligand2)]
			mcss += prot.lm.mcss.is_valid(m)

print(total, mcss)
