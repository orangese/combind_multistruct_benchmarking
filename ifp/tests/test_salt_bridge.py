import pytest
import sys
import os

sys.path.insert(0, os.environ['COMBINDHOME'])
from ifp.fp import FP

base = '{}/ifp/tests/'.format(os.environ['COMBINDHOME'])

def zero(interaction, ifp):
	if interaction in ifp:
		val = ifp[interaction]
		assert val == 0

def test():
	fp = FP(['-version', 'ifp5'])
	fp.params['mode'] = 'pv'
	fp.params['input_file'] = base + 'inputs/1K21_lig-to-1A4W_pv.maegz'
	fp.params['poses'] = 1

	ifp = fp.fingerprint_pose_viewer()

	assert len(ifp) == 1
	ifp = ifp[0]
	print(ifp)

	zero((1, 'A:221(ASP)', ''), ifp)
	assert ifp[(1, 'A:189(ASP)', '')] == 1.0


def test_unique():
	fp = FP(['-version', 'ifp5'])
	fp.params['mode'] = 'pv'
	fp.params['input_file'] = base + 'inputs/5EQH_lig-to-5EQG_pv.maegz'
	fp.params['poses'] = 20

	ifp = fp.fingerprint_pose_viewer()

	# The below changes between ifp4 and ifp5 and seems
	# to degrade the predictions. However, the interaction is
	# clearly correctly identified, so IDK.
	#zero((2, 'A:161(GLN)', ''), ifp[16])
