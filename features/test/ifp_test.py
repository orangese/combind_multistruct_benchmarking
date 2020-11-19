import pytest
import ifp
import gzip
from rdkit.Chem.rdmolfiles import MaeMolSupplier

with gzip.open('test/pv.maegz') as fp:
    mols =  MaeMolSupplier(fp, removeHs=False)
    protein = next(mols)
    ligands = list(mols)

settings = {'version'           : 'rd1',
            'level'             : 'residue',
            'hbond_dist_opt'    : 3.5,
            'hbond_dist_cut'    : 4.0,
            'hbond_angle_opt'   : 60.0,
            'hbond_angle_cut'   : 90.0,
            'sb_dist_opt'       : 4.0,
            'sb_dist_cut'       : 5.0,
            'contact_scale_opt' : 1.25,
            'contact_scale_cut' : 1.75,
            'pipi_dist_opt'     : 7.0,
            'pipi_dist_cut'     : 8.0}

settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
def test_version():
    import rdkit
    assert rdkit.__version__ == '2020.03.1'

def test_hydrogenbond():
    i = ifp.hbond_compute(protein, ligands[0], settings)
    assert len(i) == 3

def test_saltbridge_none():
    i = ifp.saltbridge_compute(protein, ligands[0], settings)
    assert len(i) == 0

def test_saltbridge_one():
    i = ifp.saltbridge_compute(protein, ligands[3], settings)
    assert len(i) == 1

def test_contact():
    ifp.contact_compute(protein, ligands[0], settings)

def test_pipi_tstack():
    i = ifp.pipi_compute(protein, ligands[0], settings)
    assert len(i) == 1
    i = ifp.pipi_compute(protein, ligands[3], settings)
    assert len(i) == 1

def test_pipi_pstack():
    i = ifp.pipi_compute(protein, ligands[173], settings)
    print(i)
    assert len(i) == 2

    i = ifp.pipi_compute(protein, ligands[180], settings)
    print(i)
    assert len(i) == 2

# with gzip.open('test/pv.maegz') as fp:
#     mols =  MaeMolSupplier(fp, removeHs=False)
#     protein = next(mols)
#     ifp.fingerprint_poseviewer(protein, mols, 100, settings)
