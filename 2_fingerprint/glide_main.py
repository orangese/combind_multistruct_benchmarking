#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
from os import listdir
from os.path import isfile, join
import itertools

def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

SCRIPT = '/share/PI/rondror/$USER/combind/2_fingerprint/fuzzyifp.py'
SCHRODINGER = '/share/PI/rondror/software/schrodinger2017-1/run'
DATA = '/scratch/PI/rondror/docking_data/'
XGLIDE = '../../xglide'

output_dir = sys.argv[1]
datasets = sys.argv[2:]

for dataset in datasets:    
    os.system('mkdir -p {}{}/ifp/{}'.format(DATA, dataset, output_dir))
    os.chdir('{}{}/ifp/{}'.format(DATA, dataset, output_dir))
    if dataset.split('/')[0] == 'pdbbind':
        title = dataset.split('/')[1]
    else:
        title = dataset
    empty = 0
    for f in listdir('.'):
        if f.split('.')[-1] == 'fp' and os.path.getsize(f) == 0: 
            empty += 1
            os.system('rm -f {}.fp {}.out'.format(f.split('.')[0], f.split('.')[0]))
    print 'deleted {} empty fingerprint files'.format(empty)
    
    glideFolders = [f for f in listdir(XGLIDE) if not os.path.exists(f+'.fp')]
    print '{}: {} of {} xglide folders have not been fingerprinted. submitting these...'.format(dataset, len(glideFolders), len(listdir(XGLIDE)))
    
    for i, folderGroup in enumerate(grouper(6, glideFolders)): #Fingerprint jobs are grouped into batches of 6
        with open('{}{}.sh'.format(title, i), 'w') as f: #Write 6 jobs to a script
            f.write("#!/bin/bash\n")
            for folder in folderGroup:
                if folder != None:
                    verbose = '{}.out'.format(folder)
                    pose_viewer = '{}/{}/{}_pv.maegz'.format(XGLIDE, folder, folder)
                    f.write('{} {} -receptor {} -input pose_viewer -verbose {}.out > {}.fp &\n'.format(SCHRODINGER, SCRIPT, pose_viewer, folder, folder))
                    #f.write(SCHRODINGER + " " + SCRIPT + " -receptor " + pose_viewer + " -input pose_viewer -verbose " + verbose_output + " > " + ifpDir+"/"+folder+".fp &\n")

            f.write("wait\n") #Wait for all forks for finish

        os.system('sbatch --time=04:00:00 -c 6 -p rondror {}{}.sh'.format(title, i))
