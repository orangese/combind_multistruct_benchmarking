import unittest
import numpy as np
from density_estimate import DensityEstimate
from statistics import merge_dicts_of_lists

class TestDensityEstimate(unittest.TestCase):
	def test_average(self):
		de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
		self.assertEqual(de1.fx.shape,  (3,))
		de1.fx = np.array([0, 0, 0])

		de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
		self.assertEqual(de2.fx.shape, (3,))
		de2.fx = np.array([1, 1, 1])

		avg = de1._average(de2)

		self.assertTrue(np.all(avg.fx == [0.75, 0.75, 0.75]), avg.fx)

	def test_average_zero(self):
		de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
		self.assertEqual(de1.fx.shape,  (3,))
		de1.fx = np.array([0, 0, 0])

		de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
		self.assertEqual(de2.fx.shape, (3,))
		de2.fx = np.array([1, 1, 1])
 
		avg = de1._average(de2)

		self.assertTrue(np.all(avg.fx == [1, 1, 1]), avg.fx)

	def test_merge(self):
		de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
		self.assertEqual(de1.fx.shape,  (3,))
		de1.fx = np.array([0, 0, 0])

		de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 15)
		self.assertEqual(de2.fx.shape, (3,))
		de2.fx = np.array([1, 1, 1])

		de3 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
		self.assertEqual(de3.fx.shape, (3,))
		de3.fx = np.array([2, 2, 2])

		merged = DensityEstimate.merge([de1, de2, de3], False)

		self.assertTrue(np.all(merged.fx == [0.75, 0.75, 0.75]), merged.fx)

	def test_merge_weighted(self):
		de1 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 5)
		self.assertEqual(de1.fx.shape,  (3,))
		de1.fx = np.array([0, 0, 0])

		de2 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 10)
		self.assertEqual(de2.fx.shape, (3,))
		de2.fx = np.array([1, 1, 1])

		de3 = DensityEstimate(points = 3, domain = (0, 2), n_samples = 0)
		self.assertEqual(de3.fx.shape, (3,))
		de3.fx = np.array([2, 2, 2])

		merged = DensityEstimate.merge([de1, de2, de3])

		self.assertTrue(np.all(merged.fx == [0.5, 0.5, 0.5]), merged.fx)

class TestMergeDicts(unittest.TestCase):
	def test_mismatched(self):
		stats1  = {}
		stats2  = {'native': {}, 'decoy': {}}

		self.assertRaises(AssertionError, merge_dicts_of_lists, stats1, stats2)

	def test_both_empty(self):
		stats1  = {'native': {}, 'decoy': {}}
		stats2  = {'native': {}, 'decoy': {}}

		merged = merge_dicts_of_lists(stats1, stats2)

		self.assertEqual(merged, {'native': {}, 'decoy': {}}, merged)

	def test_one_empty(self):
		stats1  = {'native': {'pipi': 1}, 'decoy': {}}
		stats2  = {'native': {}, 'decoy': {}}

		merged = merge_dicts_of_lists(stats1, stats2)

		self.assertEqual(merged, {'native': {'pipi': [1]}, 'decoy': {}}, merged)

	def test_full(self):
		stats1  = {'native': {'pipi': 1}, 'decoy': {}}
		stats2  = {'native': {'pipi': [2]}, 'decoy': {}}

		merged = merge_dicts_of_lists(stats1, stats2)

		self.assertEqual(merged, {'native': {'pipi': [1, 2]}, 'decoy': {}}, merged)

if __name__ == '__main__':
    unittest.main()
