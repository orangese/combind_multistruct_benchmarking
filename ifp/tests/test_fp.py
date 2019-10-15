import pytest
import sys
import os

sys.path.insert(0, os.environ['COMBINDHOME'])
import ifp.fp

def _test(name, terms):
	root = '{}/ifp/tests'.format(os.environ['COMBINDHOME'])
	args = ('ifp -poses 1 -version ifp5 -mode pv '
	        '-input_file {0:}/inputs/{1:}.maegz '
	        '-output_file {0:}/outputs/{1:}.fp').format(root, name)
	ifp.fp.FP(args.split())

	with open('{}/outputs/{}.fp'.format(root, name)) as fp:
		test = fp.read()
		for term in terms:
			assert term in test

def test_h_index_repeated():
	_test('3IPH_lig-to-1KV1_pv', ['2-A:168(ASP)-=1'])

def test_ligand_with_imine():
	_test('1F5L_lig-to-1C5X_pv', [])

def test_ligand_formal_charge_of_minus_two():
	_test('1OWH_lig-to-1C5X_pv', [])

def test_longest_reasonable_saltbridge():
	_test('4U16_lig-to-4DAJ_pv', [])

def test_fused_ring_pipi():
	_test('2HVC_lig-to-1T5Z_pv', [])

def test_fused_rings_and_disjoint_rings():
	_test('4AQC_lig-to-2B7A_pv', [])

