import pytest
from statistics import merge_dicts_of_lists

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