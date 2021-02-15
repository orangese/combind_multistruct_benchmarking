import os
import shutil

old_name = 'rd1_shape_5'
new_name = 'all_rd1_shape_5'

#os.mkdir(new_name)

for i in range(1, 5):
  new_root = '{}/{}'.format(new_name, i)
  old_root = '{}/{}'.format(old_name, i)
    
  os.mkdir(new_root)
  shutil.copy('{}/binder.smi'.format(old_root), new_root)

  shutil.copytree('{}/bpp'.format(old_root), '{}/bpp'.format(new_root))
