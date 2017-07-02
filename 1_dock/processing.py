import os
from multiprocessing import Pool
import slurm

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(fileName):
    options = "-WAIT -fillsidechains -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini"  
    command = "{}/utilities/prepwizard {} -j temp-{} ../stripped/{}.mae {}.mae".format(SCHRODINGER, options, fileName, fileName, fileName)
    slurm.salloc(command, "1", "5:00:00") #(jobString, numProcessors, timeLimit)

def process():   
    currentFiles = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
    currentFiles = map(lambda x: os.path.splitext(x)[0], currentFiles)

    #Each thread has a very lightweight job (waiting on the salloc signal) so just launch len(currentFiles) threads
    pool = Pool(len(currentFiles)) 
    pool.map(processHelper, currentFiles)

