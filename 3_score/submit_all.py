import os
import sys

sys.path.append('../notebooks')
from load_data import *

DATA = '/scratch/PI/rondror/docking_data'
pdbbind = ['pdbbind_combo/{}'.format(s) for s in os.listdir('{}/pdbbind_combo'.format(DATA))]
pdbbind.sort()

grids = load_grids()

for i in pdbbind:
    for g in grids[i]:
        print i, g
        #os.system('python score.py {} {}'.format(i, g))
        #break
    #break
#for i in ['B1AR_all', 'B2AR_all']:
 #   os.system('rm -rf {}/{}/scoring_output'.format(DATA, i))
 #   for g in os.listdir('{}/{}/grids'.format(DATA, i)):
 #       print i, g
        os.system('sbatch --job-name={}-{} score.sh {} {}'.format(g, i.split('_')[0], i, g))
        #print i, g
        #break
    #break
