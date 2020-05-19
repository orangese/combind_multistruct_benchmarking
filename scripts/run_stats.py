import os
import sys
sys.path.append(os.environ['COMBINDHOME'])
import utils
import config

data       = '/oak/stanford/groups/rondror/users/jpaggi/combind'
stats_root = '/oak/stanford/groups/rondror/users/jpaggi/temp/rd1'

for protein in os.listdir(data):
	if protein[0] == '.': continue
	cmd = './main.py --data {} statistics --stats_version rd1 {} {}'.format(data, stats_root, protein)
	log = '{}/{}.log'.format(stats_root, protein)
	os.system('sbatch -p rondror -t 01:00:00 -o {} --wrap="{}"'.format(log, cmd))
