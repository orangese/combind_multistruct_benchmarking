import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from settings import stats, proteins, paths

for prot in proteins:
        protein = Protein(prot, stats['stats41'], paths)
        ligands = protein.lm.get_xdocked_ligands(20)
        for lig in ligands:
                print(prot, protein.lm.st, lig.replace('_lig', ''))
