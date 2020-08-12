import pytest
from ifp import ifp
import gzip
from rdkit.Chem.rdmolfiles import MaeMolSupplier

with gzip.open('pv.maegz') as fp:
	mols =  MaeMolSupplier(fp, removeHs=False)
	protein = next(mols)
	ligands = list(mols)

_ligands = [ifp._mol_to_np(ligand) for ligand in ligands]
_protein = ifp._mol_to_np(protein)
protein_coords = ifp._atoms_to_coords(_protein)

settings = {'version'           : 'rd1',
			'level'             : 'residue',
			'hbond_dist_opt'    : 3.5,
			'hbond_dist_cut'    : 4.0,
			'hbond_angle_opt'   : 60.0,
			'hbond_angle_cut'   : 90.0,
			'sb_dist_opt'       : 4.0,
			'sb_dist_cut'       : 5.0,
			'contact_scale_opt' : 1.25,
			'contact_scale_cut' : 1.75}

settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
settings['overall_cut'] = max(settings['hbond_dist_cut'], settings['sb_dist_cut'],
                              settings['contact_scale_cut']*2*max(settings['nonpolar'].values()))

def test_hydrogenbond():
	protein = ifp._relevent_atoms(_protein, _ligands[0], protein_coords, settings['overall_cut'])
	i = ifp.hbond_compute(protein, _ligands[0], settings)
	assert len(i) == 3

def test_saltbridge_none():
	i = ifp.saltbridge_compute(protein.GetAtoms(), ligands[0].GetAtoms(), settings)
	assert len(i) == 0

def test_saltbridge_one():
	i = ifp.saltbridge_compute(protein.GetAtoms(), ligands[3].GetAtoms(), settings)
	assert len(i) == 1

def test_contact():
	ifp.contact_compute(protein.GetAtoms(), ligands[0].GetAtoms(), settings)

def test_fingerprint_poseviewer():
	ifp.fingerprint_poseviewer(protein, ligands, settings)
