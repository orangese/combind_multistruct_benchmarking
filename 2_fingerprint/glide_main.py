import os
import sys
from os import listdir
from os.path import isfile, join
import itertools

def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

SCRIPT = '/share/PI/rondror/docking_code/2_fingerprint/fuzzyifp.py'
SCHRODINGER = '/share/PI/rondror/software/schrodinger2017-1/run'

glideDir = sys.argv[1]
ifpDir = sys.argv[2]

os.system("mkdir -p " + ifpDir)

glideFolders = [f for f in listdir(glideDir)]

for groupNum, folderGroup in enumerate(grouper(6, glideFolders)): #Fingerprint jobs are grouped into batches of 6
    with open("fifpjob_batch_" + str(groupNum) + ".sh", "w") as f: #Write 6 jobs to a script
        f.write("#!/bin/bash\n")
        for folder in folderGroup:
            if folder != None:
                f.write(SCHRODINGER + " " + SCRIPT + " -receptor " + glideDir+"/"+folder+"/"+folder+"_pv.maegz -input pose_viewer > " + ifpDir+"/"+folder+".fp &\n")

        f.write("wait\n") #Wait for all forks for finish

    os.system("sbatch --time=50:00:00 -c 6 -p rondror " + "fifpjob_batch_" + str(groupNum) + ".sh") #Submit the script
