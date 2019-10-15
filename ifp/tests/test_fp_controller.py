import pytest
import os
import sys

from ifp.fp_controller import compute_fp, _structure_fp, _compute_fp


class MockLigandManager:
	def __init__(self, tmpdir):
		self.pdb = ['6DDF_lig']
		self.st = '5C1M'
		self.params = {'ifp_version': 'ifp5',
		               'max_poses': 1}
		self.tmpdir = tmpdir

	def path(self, name, extras=None):
		if name == 'IFP':
			return self.tmpdir+'ifp/ifp5/{ligand}-to-5C1M-confgen_es11.fp'.format(**extras)

		elif name == 'IFP_ROOT':
			return self.tmpdir+'ifp/ifp5'
		elif name == 'DOCK_PV':
			return self.tmpdir+'confgen_es11/{ligand}-to-5C1M'.format(**extras)

	def chembl(self):
		return ['CHEMBL1_lig', 'CHEMBL2_lig']

	def docked(self, ligands):
		return ligands

class TestComputeFP:
	def test_compute_fp_absent(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager(tmpdir)
		mock_exists = mocker.patch('os.path.exists')
		mock_structure_fp = mocker.patch('ifp.fp_controller._structure_fp')
		mock_compute_fp = mocker.patch('ifp.fp_controller._compute_fp')
		
		compute_fp(lm)

		mock_exists.assert_any_call(tmpdir+'ifp/ifp5/CHEMBL1_lig-to-5C1M-confgen_es11.fp')
		mock_exists.assert_any_call(tmpdir+'ifp/ifp5/CHEMBL2_lig-to-5C1M-confgen_es11.fp')
		mock_exists.assert_any_call(tmpdir+'ifp/ifp5/6DDF_lig-to-5C1M-confgen_es11.fp')
		assert mock_exists.call_count == 3

		mock_structure_fp.assert_called_with(lm)
		mock_compute_fp.assert_called_with(lm, [])

	def test_compute_fp_present(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager(tmpdir)
		mocker.patch('os.path.exists', lambda x: False)
		mock_structure_fp = mocker.patch('ifp.fp_controller._structure_fp')
		mock_compute_fp = mocker.patch('ifp.fp_controller._compute_fp')
		compute_fp(lm)

		mock_structure_fp.assert_called_with(lm)
		mock_compute_fp.assert_called_with(lm, ['6DDF_lig',
		                               	    'CHEMBL1_lig',
		                                    'CHEMBL2_lig'
		                                    ])

class TestGetFP:
	def test__compute_fp(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager(tmpdir)
		mock_system = mocker.patch('os.system')
		
		_compute_fp(lm, lm.pdb+lm.chembl())

		assert mock_system.call_count == 1
		assert os.path.exists('0fp.sh')
		with open('0fp.sh') as fp:
			lines = fp.read().split('\n')
		
		assert any(['6DDF_lig-to-5C1M-confgen_es11.fp' in line for line in lines])
		assert any(['confgen_es11/6DDF_lig-to-5C1M'    in line for line in lines])
		assert any(['CHEMBL1_lig-to-5C1M-confgen_es11.fp' in line for line in lines])
		assert any(['confgen_es11/CHEMBL1_lig-to-5C1M'    in line for line in lines])
		assert any(['CHEMBL2_lig-to-5C1M-confgen_es11.fp' in line for line in lines])
		assert any(['confgen_es11/CHEMBL2_lig-to-5C1M'    in line for line in lines])

