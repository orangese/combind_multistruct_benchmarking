"""
Tests for optimization code.

Note that these current tests are quite basic and
ought to be expanded.
"""

import pytest
import numpy as np

from score.prob_opt import PredictStructs
from score.density_estimate import DensityEstimate
from containers import Ligand, Pose


def create_ligand(name, rmsds, gscores, emodels, fps):
	ligand = Ligand(name, '', '', '')
	ligand.poses = []
	for r, g, e, f in zip(rmsds, gscores, emodels, fps):
		ligand.poses += [Pose(r, g, e, f)]
	return ligand


def create_stats():
	stats = {'native': {'sb': DensityEstimate(domain = (0, 1))},
			 'reference': {'sb': DensityEstimate(domain = (0, 1))}}
	stats['native']['sb'].fx = np.linspace(0.01, 2, stats['native']['sb'].fx.shape[0])
	stats['reference']['sb'].fx = np.linspace(2, 0.01, stats['native']['sb'].fx.shape[0])
	return stats

def basic_ps():
	ligands = {'lig1': create_ligand('lig1', [0, 1], [-2, -1.5], [0, 0], [{}, {(1, 23): 1.0}]),
			   'lig2': create_ligand('lig2', [0, 1], [-10, -3.5], [0, 0],[{(1, 23): 1.0}, {}])}
	stats = create_stats()
	features = ['sb']
	return PredictStructs(ligands, None, stats, features, 3, 1.0)

def test_get():
	ps = basic_ps()
	assert ps._get_gscore('lig1', 0) == -2
	assert ps._get_gscore('lig1', 1) == -1.5
	assert ps._get_gscore('lig2', 1) == -3.5
	assert ps._num_poses('lig1') == 2
	assert ps._num_poses('lig2') == 2

	assert ps._get_feature('sb', 'lig1', 'lig2', 0, 0) == 0
	assert ps._get_feature('sb', 'lig1', 'lig2', 0, 1) == 0
	assert ps._get_feature('sb', 'lig1', 'lig2', 1, 0) == 1.0
	assert ps._get_feature('sb', 'lig1', 'lig2', 1, 1) == 0
	assert ps._get_feature('sb', 'lig2', 'lig1', 0, 0) == 0
	assert ps._get_feature('sb', 'lig2', 'lig1', 0, 1) == 1.0
	assert ps._get_feature('sb', 'lig2', 'lig1', 1, 0) == 0
	assert ps._get_feature('sb', 'lig2', 'lig1', 1, 1) == 0

def test_like():
	ps = basic_ps()
	def like1(p1, p2):
		return ps._likelihoods_for_pair_and_single_feature('sb',{'lig1':p1, 'lig2':p2},'lig1', 'lig2')
	def like2(p1, p2):
		return ps._likelihoods_for_pair_and_single_feature('sb',{'lig1':p1, 'lig2':p2},'lig1', 'lig2')


	assert like1(0, 0) == (0.0, 0.01, 2.0)
	assert like1(0, 1) == (0.0, 0.01, 2.0)
	assert like1(1, 0) == (1.0, 2.0, 0.01)
	assert like1(1, 1) == (0.0, 0.01, 2.0)
	assert like2(0, 0) == (0.0, 0.01, 2.0)
	assert like2(0, 1) == (0.0, 0.01, 2.0)
	assert like2(1, 0) == (1.0, 2.0, 0.01)
	assert like2(1, 1) == (0.0, 0.01, 2.0)

def test_ratio():
	ps = basic_ps()
	def like(p1, p2):
		return ps._log_likelihood_ratio_pair({'lig1':p1, 'lig2':p2},'lig1', 'lig2')

	assert like(0, 0) == np.log(0.01) - np.log(2.0)
	assert like(0, 1) == np.log(0.01) - np.log(2.0)
	assert like(1, 0) == np.log(2.0) - np.log(0.01)
	assert like(1, 1) == np.log(0.01) - np.log(2.0)

def test_partial():
	ps = basic_ps()
	def like(p1, p2, lig):
		return ps._partial_log_posterior({'lig1':p1, 'lig2':p2}, lig)

	assert like(0, 0, 'lig1') == np.log(0.01) - np.log(2.0) + 2*1.0
	assert like(0, 1, 'lig1') == np.log(0.01) - np.log(2.0) + 2*1.0
	assert like(1, 0, 'lig1') == np.log(2.0) - np.log(0.01) + 1.5*1.0
	assert like(1, 1, 'lig1') == np.log(0.01) - np.log(2.0) + 1.5*1.0

	assert like(0, 0, 'lig2') == np.log(0.01) - np.log(2.0) + 10*1.0
	assert like(0, 1, 'lig2') == np.log(0.01) - np.log(2.0) + 3.5*1.0
	assert like(1, 0, 'lig2') == np.log(2.0) - np.log(0.01) + 10.0*1.0
	assert like(1, 1, 'lig2') == np.log(0.01) - np.log(2.0) + 3.5*1.0

def test_posterior():
	ps = basic_ps()
	def like(p1, p2):
		return ps.log_posterior({'lig1':p1, 'lig2':p2})

	assert like(0, 0) == np.log(0.01) - np.log(2.0) + 2*1.0 + 10*1.0
	assert like(0, 1) == np.log(0.01) - np.log(2.0) + 2*1.0 + 3.5*1.0
	assert like(1, 0) == np.log(2.0) - np.log(0.01) + 1.5*1.0 + 10.0*1.0
	assert like(1, 1) == np.log(0.01) - np.log(2.0) + 1.5*1.0 + 3.5*1.0

def test_optimize():
	ps = basic_ps()
	def like(p1, p2):
		return ps._optimize_cluster({'lig1':p1, 'lig2':p2}, 1)

	opt = {'lig1':1, 'lig2':0}
	assert like(0, 0)[1] == opt
	assert like(0, 1)[1] == opt
	assert like(1, 0)[1] == opt
	assert like(1, 1)[1] == opt

def test_max():
	ps = basic_ps()

	opt = {'lig1':1, 'lig2':0}
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt

