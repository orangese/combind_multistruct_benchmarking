import os
import sys
sys.path.append(os.environ['COMBINDHOME'])
from settings import proteins
from score.statistics import Statistics

version = sys.argv[1]
stats_root = '/oak/stanford/groups/rondror/users/jpaggi/statistics'
root = '{}/{}'.format(stats_root, version)
os.system('mkdir -p ' + root)

for protein in proteins:
      cmd = './main.py statistics {} {} {}'.format(version, protein, root)
      log = '{}/{}.log'.format(root, protein)
      os.system('sbatch -p owners -t 01:00:00 -o {} --wrap="{}"'.format(log, cmd))