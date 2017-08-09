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

dataset = sys.argv[1]
output_dir = sys.argv[2]

if dataset == 'all':
    datasets = ['AR','MAP4K4','B2AR','CDK2','CHK1','HSP90','LPXC','TRMD']
else:
    datasets = [dataset]

for dataset in datasets:    
    glideDir = DATA + dataset + '/xglide/'
    ifpDir = DATA + dataset + '/ifp/' + output_dir + '/'

    os.system("mkdir -p " + ifpDir)
    os.chdir(ifpDir)

    empty = 0
    for f in listdir('.'):
        name, ext = f.split('.')
        if ext == 'fp':
            if os.path.getsize(f) == 0: 
                empty += 1
                os.system('rm -f {}.fp {}.out'.format(name, name))
    print 'deleted {} empty fingerprint files'.format(empty)
    
    glideFolders = [f for f in listdir(glideDir) if not os.path.exists(ifpDir+f+'.fp')]
    print '{}: {} of {} xglide folders have not been fingerprinted. submitting these...'.format(dataset, len(glideFolders), len(listdir(glideDir)))
    
    for i, folderGroup in enumerate(grouper(6, glideFolders)): #Fingerprint jobs are grouped into batches of 6
        with open(dataset + str(i) + ".sh", "w") as f: #Write 6 jobs to a script
            f.write("#!/bin/bash\n")
            for folder in folderGroup:
                if folder != None:
                    verbose_output = '{}/{}.out'.format(ifpDir, folder)
                    pose_viewer = '{}/{}/{}_pv.maegz'.format(glideDir, folder, folder)
                    f.write(SCHRODINGER + " " + SCRIPT + " -receptor " + pose_viewer + " -input pose_viewer -verbose " + verbose_output + " > " + ifpDir+"/"+folder+".fp &\n")

            f.write("wait\n") #Wait for all forks for finish

        os.system("sbatch --time=03:00:00 -c 6 -p rondror " + dataset + str(i) + ".sh") #Submit the script
