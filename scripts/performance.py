import pandas as pd
import sys

data = []
for path in sys.argv[1:]:
	_data = pd.read_csv(path)
	_data = _data.loc[_data['glide_rmsd'] != 'None']
	_data = _data.loc[_data['combind_rmsd'] != 'best=0']
	if not len(_data):
		print(path)
	data += [_data]

data = pd.concat(data)

data = data[['combind_rmsd', 'glide_rmsd', 'best_rmsd']].astype(float)

data['total'] = 1

data = data <= 2.05

print(data.sum(axis=0))
