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

glide = '/xglide/'

if dataset == 'all':
    datasets = ['AR','MAP4K4','B2AR','CDK2','CHK1','HSP90','LPXC','TRMD']
else:
    datasets = [dataset]

for dataset in datasets:
    print dataset
    
    glideDir = DATA + dataset + glide
    ifpDir = DATA + dataset + '/ifp/' + output_dir + '/'

    if os.path.exists(ifpDir):
        print dataset + ' already has that output directory. continuing...'
        continue

    os.system("mkdir " + ifpDir)
    os.chdir(ifpDir)

    glideFolders = [f for f in listdir(glideDir)]

    for groupNum, folderGroup in enumerate(grouper(6, glideFolders)): #Fingerprint jobs are grouped into batches of 6
        with open(dataset + str(groupNum) + ".sh", "w") as f: #Write 6 jobs to a script
            f.write("#!/bin/bash\n")
            for folder in folderGroup:
                if folder != None:
                    verbose_output = ifpDir + '/' + folder + '.out'
                    f.write(SCHRODINGER + " " + SCRIPT + " -receptor " + glideDir+"/"+folder+"/"+folder+"_pv.maegz -input pose_viewer -verbose " + verbose_output + " > " + ifpDir+"/"+folder+".fp &\n")

            f.write("wait\n") #Wait for all forks for finish

        os.system("sbatch --time=01:00:00 -c 6 -p rondror " + dataset + str(groupNum) + ".sh") #Submit the script
