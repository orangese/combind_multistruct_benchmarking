import os
import sys
sys.path.append(os.environ['COMBINDHOME'])
import utils
import config

data       = '/oak/stanford/groups/rondror/users/jpaggi/combind'
stats_root = '/oak/stanford/groups/rondror/users/jpaggi/temp/rd1'

# 5HT2B B1AR BACE1 CDK2 DHFR ERA F2 HSP90AA1 NR3C1 P00760 PLAU PYGM SIGMAR1 SMO AR B2AR BRD4 DAT ELANE F10 GLUT1 MGLUR5 NR3C2 PDE10A PTPN1 SLC6A4 VDR
stats_prots = ['5HT2B', 'B1AR', 'BACE1',  'CDK2', 'DHFR', 'ERA', 'F2', 'HSP90AA1',
			   'NR3C1', 'P00760', 'PLAU', 'PYGM', 'SIGMAR1', 'SMO', 'AR', 'B2AR',
			   'BRD4', 'DAT', 'ELANE', 'F10', 'GLUT1', 'MGLUR5', 'NR3C2', 'PDE10A',
			   'PTPN1', 'SLC6A4', 'VDR']

for protein in os.listdir(data):
	if protein[0] == '.': continue
	if protein not in stats_prots: continue
	cmd = './main.py --data {} statistics --stats_version rd1 {} {}'.format(data, stats_root, protein)
	log = '{}/{}.log'.format(stats_root, protein)
	os.system('sbatch -p rondror -t 01:00:00 -o {} --wrap="{}"'.format(log, cmd))
