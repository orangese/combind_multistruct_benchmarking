"""
python scripts/pick_helpers.py /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/structures/pdb.csv /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/chembl/CHEMBL*.csv /oak/stanford/groups/rondror/users/jpaggi/combind/P19491/chembl --criteria affinity
"""

import subprocess
from glob import glob

path = '/oak/stanford/groups/rondror/users/jpaggi/combind/{}/chembl'
chembl = 'stats_data/systems.txt'

cmd = 'sbatch -p rondror -t 02:00:00 --wrap="python ~/combind/scripts/chembl.py {}"'

with open(chembl) as fp:
	for line in fp:
		protein, chembl_id = line.strip().split()

		if glob(path.format(protein) + '/CHEMBL*'):
			continue

		cwd = path.format(protein)
		

		if protein == 'Q05586-Q12879':
			_cmd = cmd.format(chembl_id + ' --protein-complex')

		else:
			_cmd = cmd.format(chembl_id)

		print(cwd, _cmd)
		subprocess.run(_cmd, shell=True, cwd=cwd)


cmd = 'python scripts/pick_helpers.py {0}/structures/pdb.csv {0}/chembl/CHEMBL*.csv {0}/chembl --criteria {1}'

with open(chembl) as fp:
	for line in fp:
		protein, chembl_id = line.strip().split()

		if glob(path.format(protein) + '/CHEMBL*'):
			continue

		root = '{}/{}'.format(root, protein)
		_cmd = cmd.format(root, )

		print(cwd, _cmd)
		subprocess.run(_cmd, shell=True)


