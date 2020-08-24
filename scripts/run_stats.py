"""
python scripts/run_stats.py rd1 /oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems /oak/stanford/groups/rondror/users/jpaggi/temp/shape
"""

import os
import sys
sys.path.append(os.environ['COMBINDHOME'])
import utils
import config

version    = sys.argv[1]
data       = sys.argv[2]
stats_root = sys.argv[3]

CMD = './main.py --data {} statistics --stats-version {} {} {}'
for protein in os.listdir(data):
	if protein[0] == '.': continue
	cmd = CMD.format(data, version, stats_root, protein)
	log = '{}/{}.log'.format(stats_root, protein)
	os.system('sbatch -p rondror -t 01:00:00 -o {} --wrap="{}"'.format(log, cmd))
