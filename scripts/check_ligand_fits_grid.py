import sys
from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid


with StructureReader(sys.argv[1]) as st:
	ligand = list(st)[0]

with open(sys.argv[2]) as fp:
	grid = {line.strip().split()[0]: line.strip().split()[1]
			for line in fp}
center = [float(c) for c in grid['GRID_CENTER'].split(',')]
inner  = [float(c)/2 for c in grid['INNERBOX'].split(',')]
outer  = [float(c)/2 for c in grid['OUTERBOX'].split(',')]

for atom in ligand.atom:
	x, y, z = atom.xyz

	d = abs(x - center[0]) - outer[0]
	if d >= 0:
		print('OUTERBOX exceeded in x-dimension by {} for {}.'.format(d, sys.argv[1]))

	d = abs(y - center[1]) - outer[1]
	if d >= 0:
		print('OUTERBOX exceeded in y-dimension by {} for {}.'.format(d, sys.argv[1]))

	d = abs(z - center[2]) - outer[2]
	if d >= 0:
		print('OUTERBOX exceeded in z-dimension by {} for {}.'.format(d, sys.argv[1]))

x, y, z, _ = get_centroid(ligand)

d = abs(x - center[0]) - inner[0]
if d >= 0:
	print('INNERBOX exceeded in x-dimension by {} for {}.'.format(d, sys.argv[1]))

d = abs(y - center[1]) - inner[1]
if d >= 0:
	print('INNERBOX exceeded in y-dimension by {} for {}.'.format(d, sys.argv[1]))

d = abs(z - center[2]) - inner[2]
if d >= 0:
	print('INNERBOX exceeded in z-dimension by {} for {}.'.format(d, sys.argv[1]))

