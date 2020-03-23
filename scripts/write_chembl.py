import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
from settings import stats, proteins, paths

for prot in proteins:
        protein = Protein(prot, stats['stats41'], paths)
        ligands = protein.lm.get_xdocked_ligands(20)
        for lig in ligands:
                helpers = protein.lm.load_helpers()
                chembl = set(__v
                             for v in helpers.values()
                             for _v in v.values()
                             for __v in _v)
                for chembl_lig in chembl:
                        print(prot, protein.lm.st, chembl_lig.replace('_lig', ''))
