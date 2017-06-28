import os
import multiprocessing as mp

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
def process(pool):
    currentFiles = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
    currentFiles = map(lambda x: os.path.splitext(x)[0], currentFiles)
    pool.map(processHelper, currentFiles)

def processHelper(fileName):
    os.system(SCHRODINGER + "/utilities/prepwizard -WAIT -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini -j temp-" + fileName +" ../stripped/"+fileName+".mae "+fileName+".mae")
