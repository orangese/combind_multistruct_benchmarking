import sys
from schrodinger.structure import StructureReader


pv = sys.argv[1]

rank = -1

with StructureReader(pv) as st:
	for pose in st:
		if 'i_glide_XP_PoseRank' not in pose.property:
			assert rank == -1
			continue
		assert pose.property['i_glide_XP_PoseRank'] > rank
		rank = pose.property['i_glide_XP_PoseRank']
		print(rank)