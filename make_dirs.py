from shared_paths import shared_paths
import os

dirs = ['A2AR',
 'AR',
 'B2AR',
 'BRAF',
 'CDK2',
 'ERA',
 'GCR',
 'JAK2',
 # 'MR',
 # 'P00734',
 # 'P07900',
 # 'P18031',
 'PLK1'
 ]

for system in dirs:
	original_path = "{}/{}".format(shared_paths['read_data'], system)
	assert os.path.exists(original_path), "{} doesn't exist in read data".format(system)
	new_path = "{}/{}/dude/".format(shared_paths['write_data'], system)
	os.system("mkdir -p {}".format(new_path))
