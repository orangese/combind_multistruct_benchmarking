import os
import shutil

old_name = 'rd1_all_5'
new_name = 'rd1_shape_5'

#os.mkdir(new_name)

for i in range(0, 1):
	new_root = '{}/{}'.format(new_name, i)
	old_root = '{}/{}'.format(old_name, i)
	
	os.mkdir(new_root)
	shutil.copy('{}/binder.smi'.format(old_root), new_root)
	shutil.copy('{}/similarity.csv'.format(old_root), new_root)
	shutil.copy('{}/shape.csv'.format(old_root), new_root)
	shutil.copytree('{}/shape'.format(old_root), '{}/shape'.format(new_root))

	os.mkdir('{}/bpp'.format(new_root))
	shutil.copytree('{}/bpp/docking'.format(old_root), '{}/bpp/docking'.format(new_root))
	shutil.copytree('{}/bpp/ligands'.format(old_root), '{}/bpp/ligands'.format(new_root))
	# shutil.copytree('{}/bpp/shape'.format(old_root), '{}/bpp/shape'.format(new_root))
	# shutil.copytree('{}/bpp/ifp-pair'.format(old_root), '{}/bpp/ifp-pair'.format(new_root))
	# shutil.copytree('{}/bpp/mcss'.format(old_root), '{}/bpp/mcss'.format(new_root))
