"""
Tests for the LigPair class.
"""

import pytest
from score.pairs import LigPair
from containers import Ligand, Pose

def create_ligand(name, fps):
	ligand = Ligand(name, '', '', '')
	ligand.poses = []
	for fp in fps:
		ligand.poses += [Pose(0, 0, 0, fp)]
	return ligand

def test_empty():
	lig1 = create_ligand('lig1', [{}]*5)
	lig2 = create_ligand('lig2', [{}]*5)

	lp = LigPair(lig1, lig2, ['sb', 'hbond', 'mcss'], None, 4)

	assert len(lp.pose_pairs) == 16
	assert lp.feat_map == {'sb': (0.0, 0.0),
						   'hbond': (0.0, 0.0),
						   'mcss': (float('inf'), -float('inf'))}

	for i in range(4):
		for j in range(4):
			assert lp.get_feature('sb', i, j) is None
			assert lp.get_feature('hbond', i, j) is None
			assert lp.get_feature('mcss', i, j) is None
			with pytest.raises(KeyError):
				lp.get_feature('contact', i, j)

def test_max_one():
	lig1 = create_ligand('lig1', [{(1, 23): 1.0}, {(1, 23): 0.5}, {}])
	lig2 = create_ligand('lig2', [{}, {(1, 23): 0.5}, {(1, 23): 1.0}])

	lp = LigPair(lig1, lig2, ['sb', 'hbond', 'mcss'], None, 4)

	assert len(lp.pose_pairs) == 9
	assert lp.feat_map == {'sb': (0.0, 1.0),
						   'hbond': (0.0, 0.0),
						   'mcss': (float('inf'), -float('inf'))}

	assert lp.get_feature('sb', 0, 0) == 0.0
	assert lp.get_feature('sb', 1, 0) == 0.0
	assert lp.get_feature('sb', 2, 0) == 0.0
	assert lp.get_feature('sb', 0, 1) == 0.5**0.5
	assert lp.get_feature('sb', 1, 1) == (0.5*0.5)**0.5
	assert lp.get_feature('sb', 2, 1) == 0.0
	assert lp.get_feature('sb', 0, 2) == 1.0
	assert lp.get_feature('sb', 1, 2) == 0.5**0.5
	assert lp.get_feature('sb', 2, 2) == 0.0

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) is None
			assert lp.get_feature('mcss', i, j) is None
			with pytest.raises(KeyError):
				lp.get_feature('contact', i, j)

def test_max_less_than_one():
	lig1 = create_ligand('lig1', [{(1, 23): 0.9}, {(1, 23): 0.5}, {}])
	lig2 = create_ligand('lig2', [{}, {(1, 23): 0.5}, {(1, 23): 1.0}])

	lp = LigPair(lig1, lig2, ['sb', 'hbond', 'mcss'], None, 4)

	assert len(lp.pose_pairs) == 9
	assert lp.feat_map == {'sb': (0.0, 0.9**0.5),
						   'hbond': (0.0, 0.0),
						   'mcss': (float('inf'), -float('inf'))}

	assert lp.get_feature('sb', 0, 0) == 0.0
	assert lp.get_feature('sb', 1, 0) == 0.0
	assert lp.get_feature('sb', 2, 0) == 0.0
	assert lp.get_feature('sb', 0, 1) == (0.9*0.5)**0.5
	assert lp.get_feature('sb', 1, 1) == (0.5*0.5)**0.5
	assert lp.get_feature('sb', 2, 1) == 0.0
	assert lp.get_feature('sb', 0, 2) == 0.9**0.5
	assert lp.get_feature('sb', 1, 2) == 0.5**0.5
	assert lp.get_feature('sb', 2, 2) == 0.0

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) is None
			assert lp.get_feature('mcss', i, j) is None
			with pytest.raises(KeyError):
				lp.get_feature('contact', i, j)

def test_max_greater_than_one():
	lig1 = create_ligand('lig1', [{(1, 23): 1.0, (1, 20): 1.0}, {(1, 23): 0.5}, {(1, 20): 1.0}])
	lig2 = create_ligand('lig2', [{}, {(1, 23): 0.5}, {(1, 23): 1.0, (1, 20): 1.0}])

	lp = LigPair(lig1, lig2, ['sb', 'hbond', 'mcss'], None, 4)

	assert len(lp.pose_pairs) == 9
	assert lp.feat_map == {'sb': (0.0, 2.0),
						   'hbond': (0.0, 0.0),
						   'mcss': (float('inf'), -float('inf'))}

	assert lp.get_feature('sb', 0, 0) == 0.0
	assert lp.get_feature('sb', 1, 0) == 0.0
	assert lp.get_feature('sb', 2, 0) == 0.0
	assert lp.get_feature('sb', 0, 1) == 0.5**0.5 / 2.0
	assert lp.get_feature('sb', 1, 1) == (0.5*0.5)**0.5 / 2.0
	assert lp.get_feature('sb', 2, 1) == 0.0
	assert lp.get_feature('sb', 0, 2) == 1.0
	assert lp.get_feature('sb', 1, 2) == 0.5**0.5 / 2.0
	assert lp.get_feature('sb', 2, 2) == 0.5

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) is None
			assert lp.get_feature('mcss', i, j) is None
			with pytest.raises(KeyError):
				lp.get_feature('contact', i, j)