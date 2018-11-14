import pytest
import os
from score_controller import score

def test_standard_pdb(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')

	score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], False,
	      2.5, ['mcss', 'pipi'], False, None)

	assert os.path.exists('2.5-mcss_pipi')

	assert os.path.exists('2.5-mcss_pipi/settings.py')
	with open('2.5-mcss_pipi/settings.py') as fp:
		txt = fp.read()
		assert 'alpha=2.5' in txt
		assert 'use_crystal_pose=False' in txt
		assert "k_list=['mcss', 'pipi']" in txt
		assert "chembl=False" in txt
		assert "num_poses=100" in txt

	assert os.path.exists('2.5-mcss_pipi/run.sh')
	with open('2.5-mcss_pipi/run.sh') as fp:
		txt = fp.read()
		assert 'stats 5C1M MOR 4DKL_lig 6DDF_lig' in txt

def test_crystal_pdb(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')

	score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], True,
	      2.5, ['mcss', 'pipi'], False, None)

	assert os.path.exists('2.5-mcss_pipi')

	assert os.path.exists('2.5-mcss_pipi/settings.py')
	with open('2.5-mcss_pipi/settings.py') as fp:
		txt = fp.read()
		assert 'alpha=5.0' in txt
		assert 'use_crystal_pose=True' in txt
		assert "k_list=['mcss', 'pipi']" in txt
		assert "chembl=False" in txt
		assert "num_poses=100" in txt

	assert os.path.exists('2.5-mcss_pipi/run.sh')
	with open('2.5-mcss_pipi/run.sh') as fp:
		txt = fp.read()
		assert 'stats 5C1M MOR 4DKL_lig 6DDF_lig' in txt

def test_crystal_only_pdb(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')

	score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], True,
	      2.5, ['mcss', 'pipi'], True, None)

	assert os.path.exists('2.5-mcss_pipi')

	assert os.path.exists('2.5-mcss_pipi/settings.py')
	with open('2.5-mcss_pipi/settings.py') as fp:
		txt = fp.read()
		assert 'alpha=2.5' in txt
		assert 'use_crystal_pose=True' in txt
		assert "k_list=['mcss', 'pipi']" in txt
		assert "chembl=False" in txt
		assert "num_poses=100" in txt

	assert os.path.exists('2.5-mcss_pipi/run.sh')
	with open('2.5-mcss_pipi/run.sh') as fp:
		txt = fp.read()
		assert 'stats 5C1M MOR 6DDF_lig' in txt
		assert 'stats 5C1M MOR 4DKL_lig' in txt


def test_standard_chembl(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')

	score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], False,
	      2.5, ['mcss', 'pipi'], False, ('best_affinity', 5))

	assert os.path.exists('5-2.5-mcss_pipi')

	assert os.path.exists('5-2.5-mcss_pipi/settings.py')
	with open('5-2.5-mcss_pipi/settings.py') as fp:
		txt = fp.read()
		assert 'alpha=12.5' in txt
		assert 'use_crystal_pose=False' in txt
		assert "k_list=['mcss', 'pipi']" in txt
		assert "chembl=True" in txt
		assert "num_poses=100" in txt
		assert 'chembl_file="best_affinity.txt"' in txt

	assert os.path.exists('5-2.5-mcss_pipi/run.sh')
	with open('5-2.5-mcss_pipi/run.sh') as fp:
		txt = fp.read()
		assert 'stats 5C1M MOR 6DDF_lig' in txt
		assert 'stats 5C1M MOR 4DKL_lig' in txt

def test_crystal_chembl(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')

	score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], True,
	      2.5, ['mcss', 'pipi'], False, ('best_affinity', 5))

	assert os.path.exists('5-2.5-mcss_pipi')

	assert os.path.exists('5-2.5-mcss_pipi/settings.py')
	with open('5-2.5-mcss_pipi/settings.py') as fp:
		txt = fp.read()
		assert 'alpha=15.0' in txt
		assert 'use_crystal_pose=True' in txt
		assert "k_list=['mcss', 'pipi']" in txt
		assert "chembl=True" in txt
		assert "num_poses=100" in txt
		assert 'chembl_file="best_affinity.txt"' in txt

	assert os.path.exists('5-2.5-mcss_pipi/run.sh')
	with open('5-2.5-mcss_pipi/run.sh') as fp:
		txt = fp.read()
		assert 'stats 5C1M MOR 6DDF_lig' in txt
		assert 'stats 5C1M MOR 4DKL_lig' in txt

def test_crystal_only_chembl(tmpdir, mocker):
	os.chdir(tmpdir)
	mocker.patch('os.system')
	with pytest.raises(AssertionError) as excinfo:
		score('stats', '5C1M', 'MOR', ['4DKL_lig', '6DDF_lig'], True,
	          2.5, ['mcss', 'pipi'], True, ('best_affinity', 5))
	assert 'AssertionError' in str(excinfo), excinfo
