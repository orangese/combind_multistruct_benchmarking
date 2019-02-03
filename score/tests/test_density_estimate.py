import pytest
import os
import sys
import numpy as np

sys.path.insert(0, os.environ['COMBINDHOME'])
from score.density_estimate import DensityEstimate
from score.statistics import merge_dicts_of_lists

def test_average():
	de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
	assert de1.fx.shape == (3,)
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
	assert de2.fx.shape == (3,)
	de2.fx = np.array([1, 1, 1])

	avg = de1._average(de2)

	assert np.all(avg.fx == [0.75, 0.75, 0.75])

def test_average_zero():
	de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
	assert de1.fx.shape == (3,)
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
	assert de2.fx.shape == (3,)
	de2.fx = np.array([1, 1, 1])

	avg = de1._average(de2)

	assert np.all(avg.fx == [1, 1, 1])

def test_merge():
	de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
	assert de1.fx.shape == (3,)
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
	assert de2.fx.shape == (3,)
	de2.fx = np.array([1, 1, 1])

	de3 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
	assert de3.fx.shape == (3,)
	de3.fx = np.array([2, 2, 2])

	merged = DensityEstimate.merge([de1, de2, de3], False)

	assert np.all(merged.fx == [0.75, 0.75, 0.75])

def test_merge_weighted():
	de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
	assert de1.fx.shape == (3,)
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 10)
	assert de2.fx.shape == (3,)
	de2.fx = np.array([1, 1, 1])

	de3 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
	assert de3.fx.shape == (3,)
	de3.fx = np.array([2, 2, 2])

	merged = DensityEstimate.merge([de1, de2, de3])

	assert np.all(merged.fx == [0.5, 0.5, 0.5])
