import pytest
from shared_paths import shared_paths
from fp_controller import compute_fp, structure_fp, get_fp
import os

class MockLigandManager:
	def __init__(self):
		self.pdb = ['6DDF_lig']
		self.st = '5C1M'

	def chembl(self):
		return ['CHEMBL1_lig', 'CHEMBL2_lig']

	def docked(self, ligands):
		return ligands

class TestComputeFP:
	def test_compute_fp_absent(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager()
		mocker.patch.dict(shared_paths, {'docking': 'confgen_es11'})
		mock_exists = mocker.patch('os.path.exists')
		mock_structure_fp = mocker.patch('fp_controller.structure_fp')
		mock_get_fp = mocker.patch('fp_controller.get_fp')
		
		compute_fp(lm)

		mock_exists.assert_any_call('ifp/ifp4/CHEMBL1_lig-to-5C1M-confgen_es11.fp')
		mock_exists.assert_any_call('ifp/ifp4/CHEMBL2_lig-to-5C1M-confgen_es11.fp')
		mock_exists.assert_any_call('ifp/ifp4/6DDF_lig-to-5C1M-confgen_es11.fp')
		assert mock_exists.call_count == 3

		mock_structure_fp.assert_called_with(lm)
		mock_get_fp.assert_called_with(lm, [], False)

	def test_compute_fp_present(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager()
		mocker.patch('os.path.exists', lambda x: False)
		mock_structure_fp = mocker.patch('fp_controller.structure_fp')
		mock_get_fp = mocker.patch('fp_controller.get_fp')
		compute_fp(lm)

		mock_structure_fp.assert_called_with(lm)
		mock_get_fp.assert_called_with(lm, ['6DDF_lig-to-5C1M',
		                               	    'CHEMBL1_lig-to-5C1M',
		                                    'CHEMBL2_lig-to-5C1M'
		                                    ], False)

	def test_compute_fp_raw(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager()
		mocker.patch('os.path.exists', lambda x: False)
		mock_structure_fp = mocker.patch('fp_controller.structure_fp')
		mock_get_fp = mocker.patch('fp_controller.get_fp')
		compute_fp(lm, True)

		mock_structure_fp.assert_not_called()
		mock_get_fp.assert_called_with(lm, ['6DDF_lig-to-5C1M'], True)


class TestGetFP:
	def test_get_fp(self, tmpdir, mocker):
		os.chdir(tmpdir)
		lm = MockLigandManager()
		mock_system = mocker.patch('os.system')
		mocker.patch.dict(shared_paths, {'docking': 'confgen_es11'})
		
		get_fp(lm, ["{}-to-{}".format(lig, lm.st)
		            for lig in lm.pdb+lm.chembl()], False)

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

