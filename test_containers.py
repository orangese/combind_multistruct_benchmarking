import pytest
from containers import Pose, Ligand, LigandManager, Protein
from settings import stats, paths


def test_protein():
    params = stats['stats21']
    protein = Protein('B1AR', params, paths)
    protein.load_docking(['2Y00_lig', '2Y03_lig'], load_fp=True, load_mcss=True)

    assert len(protein.docking) == 1
    assert len(protein.docking[protein.lm.st]) == 2