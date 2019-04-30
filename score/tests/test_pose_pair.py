import pytest
from score.pairs import PosePair
from containers import Pose
from shared_paths import shared_paths

def pose(rmsd=0.0, gscore=0.0, emodel=0.0, fp={}):
	return Pose(rmsd, gscore, emodel, fp)

def test_correct_both():
	pose1 = pose(rmsd=1.0)
	pose2 = pose(rmsd=1.4)
	pp = PosePair(pose1, pose2, 0.0)
	
	assert pp.correct() == 1.0

def test_correct_one_1():
	pose1 = pose(rmsd=1.0)
	pose2 = pose(rmsd=shared_paths['stats']['native_thresh']+1.0)
	pp = PosePair(pose1, pose2, 0.0)
	
	assert pp.correct() == 0.0

def test_correct_one_2():
	pose1 = pose(rmsd=shared_paths['stats']['native_thresh']+1.0)
	pose2 = pose(rmsd=1.4)
	pp = PosePair(pose1, pose2, 0.0)
	
	assert pp.correct() == 0.0

def test_correct_one_neither():
	pose1 = pose(rmsd=shared_paths['stats']['native_thresh']+1.0)
	pose2 = pose(rmsd=shared_paths['stats']['native_thresh']+0.1)
	pp = PosePair(pose1, pose2, 0.0)
	
	assert pp.correct() == 0.0

def test_get_feature_empty():
	pose1 = pose(fp={})
	pose2 = pose(fp={})
	pp = PosePair(pose1, pose2, 0.0)

	assert pp.get_feature('sb') == 0.0
	assert pp.get_feature('hbond') == 0.0
	assert pp.get_feature('contact') == 0.0
	assert pp.get_feature('mcss') == 0.0

def test_get_feature_single():
	pose1 = pose(fp={(1, 23): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 4.0)

	assert pp.get_feature('sb') == 1.0
	assert pp.get_feature('hbond') == 0.0
	assert pp.get_feature('contact') == 0.0
	assert pp.get_feature('mcss') == 4.0

def test_get_feature_mismatch():
	pose1 = pose(fp={(1, 2): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 1.0)

	assert pp.get_feature('sb') == 0.0
	assert pp.get_feature('hbond') == 0.0
	assert pp.get_feature('contact') == 0.0
	assert pp.get_feature('mcss') == 1.0

def test_get_feature_multiple_of_same_type():
	pose1 = pose(fp={(1, 2): 1.0, (1, 23): 1.0})
	pose2 = pose(fp={(1, 2): 0.0, (1, 23): 1.0})
	pp = PosePair(pose1, pose2, 0.0)

	assert pp.get_feature('sb') == 1.0
	assert pp.get_feature('hbond') == 0.0
	assert pp.get_feature('contact') == 0.0
	assert pp.get_feature('mcss') == 0.0
