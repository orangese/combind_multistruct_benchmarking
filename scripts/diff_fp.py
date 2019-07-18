from containers import Ligand
import sys

protein, ligand, docking, struct, ifp1, ifp2 = sys.argv[1:]

root = '/scratch/PI/rondror/combind/bpp_data/{}'.format(protein)

dock_dir = '{}/docking/{}'.format(root, docking)
ifp1_dir = '{}/ifp/{}'.format(root, ifp1)
ifp2_dir = '{}/ifp/{}'.format(root, ifp2)

lig1 = Ligand(ligand, dock_dir, ifp1_dir, struct)
lig2 = Ligand(ligand, dock_dir, ifp2_dir, struct)

lig1.load_poses(True)
lig2.load_poses(True)

for pose in range(len(lig1.poses)):
	fp1, fp2 = lig1.poses[pose].fp, lig2.poses[pose].fp

	for interaction in set(fp1.keys()).union(fp2.keys()):
		present1 = interaction in fp1 and fp1[interaction] != 0
		present2 = interaction in fp2 and fp2[interaction] != 0

		if present1 and not present2:
			print('Pose {}: {}:{} present only in 1'.format(pose, interaction, fp1[interaction]))
		elif present2 and not present1:
			print('Pose {}: {}:{} present only in 2'.format(pose, interaction, fp2[interaction]))

		if present1 and present2:
			if fp1[interaction] != fp2[interaction]:
				print('Pose {}: {}: {} in 1 and {} in 2'.format(fp1[interaction], fp2[interaction]))