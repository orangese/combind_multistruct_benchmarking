import pytest
import sys
sys.path.append('../../')
from mcss.mcss import MCSS
from schrodinger.structure import StructureReader

base = '/home/users/jpaggi/combind/mcss'

def test_compute():
	mcss = MCSS.from_string('3PPK_lig,3PRF_lig,30,24,17,0,'
	                        'c1cncc(c12)occ2Nc(c)ccccO,'
	                        'c1cncc(c12)occ2Nc(c)ccccO,False')

	init_file = '{}/tests/outputs/temp'.format(base)
	mcss_types_file = '{}/custom_types/mcss16.typ'.format(base)
	rmsd_file = '{}/tests/outputs/3PPK_lig-3PRF_lig.csv'.format(base)
	poseviewer_paths = {'3PPK_lig': '{}/tests/inputs/3PPK_lig-to-1UWH_pv.maegz'.format(base),
					    '3PRF_lig': '{}/tests/inputs/3PRF_lig-to-1UWH_pv.maegz'.format(base)}
	pv1 = list(StructureReader(poseviewer_paths[mcss.l1])
                   )[poseviewer_paths[mcss.l1][-8:] == 'pv.maegz':]
	pv2 = list(StructureReader(poseviewer_paths[mcss.l2])
	           )[poseviewer_paths[mcss.l2][-8:] == 'pv.maegz':]
	l1_idx, l2_idx = mcss._get_atom_idxss(pv1[0], pv2[0], init_file, mcss_types_file)
	
	print(mcss.l1)
	for idxs in l1_idx:
		for idx in idxs:
			print('+'.join(map(str, idx)))

	print(mcss.l2)
	for idxs in l2_idx:
		for idx in idxs:
			print('+'.join(map(str, idx)))


	mcss.write_rmsds(poseviewer_paths, init_file, mcss_types_file,
	                 rmsd_file, 2)
	with open(rmsd_file) as fp:
		print(fp.readline())
	assert False



