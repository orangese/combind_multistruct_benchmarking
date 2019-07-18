import pytest
import numpy as np
import os
import sys

from score.density_estimate import DensityEstimate
from score.statistics import merge_dicts_of_lists, merge_stats

def test_mismatched():
	stats1  = {}
	stats2  = {'native': {}, 'decoy': {}}
	with pytest.raises(AssertionError) as excinfo:
		merge_dicts_of_lists(stats1, stats2)
	assert 'AssertionError' in str(excinfo), excinfo

def test_both_empty():
	stats1  = {'native': {}, 'decoy': {}}
	stats2  = {'native': {}, 'decoy': {}}

	merged = merge_dicts_of_lists(stats1, stats2)

	assert merged == {'native': {}, 'decoy': {}}

def test_one_empty():
	stats1  = {'native': {'pipi': 1}, 'decoy': {}}
	stats2  = {'native': {}, 'decoy': {}}

	merged = merge_dicts_of_lists(stats1, stats2)

	assert merged == {'native': {'pipi': [1]}, 'decoy': {}}, merged

def test_full():
	stats1  = {'native': {'pipi': 1}, 'decoy': {}}
	stats2  = {'native': {'pipi': [2]}, 'decoy': {}}

	merged = merge_dicts_of_lists(stats1, stats2)

	assert merged == {'native': {'pipi': [1, 2]}, 'decoy': {}}

def test_merge_stats():
	de1 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 1)
	assert de1.fx.shape == (2,)
	de1.fx = np.array([0, 0])

	de2 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 2)
	assert de2.fx.shape == (2,)
	de2.fx = np.array([1, 1])

	de3 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 5)
	assert de3.fx.shape == (2,)
	de3.fx = np.array([2, 2])

	stats = {'native': {'hbond': [de1, de2, de3]}}
	merged = merge_stats(stats, False)['native']['hbond']
	assert np.all(merged.fx == np.array([1.5, 1.5]))

	stats = {'native': {'hbond': [de1, de2, de3]}}
	merged = merge_stats(stats, True)['native']['hbond']
	assert np.all(merged.fx == np.array([1, 1]))


def test_merge_stats_zero():
	de1 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 0)
	assert de1.fx.shape == (2,)
	de1.fx = np.array([0, 0])

	de2 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 2)
	assert de2.fx.shape == (2,)
	de2.fx = np.array([1, 1])

	de3 = DensityEstimate(points = 2, domain = (0, 1), n_samples = 6)
	assert de3.fx.shape == (2,)
	de3.fx = np.array([2, 2])

	stats = {'native': {'hbond': [de1, de2, de3]}}
	merged = merge_stats(stats, False)['native']['hbond']
	assert np.all(merged.fx == np.array([1.75, 1.75]))

	stats = {'native': {'hbond': [de1, de2, de3]}}
	merged = merge_stats(stats, True)['native']['hbond']
	assert np.all(merged.fx == np.array([1.5, 1.5]))
