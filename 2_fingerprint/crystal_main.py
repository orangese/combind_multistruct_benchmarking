#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
import time

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1/run"
FUZZY_SCRIPT = "/share/PI/rondror/rondror/docking_code/2_make_fingerprints/fuzzyifp.py"

DATA = "/scratch/PI/rondror/docking_data/"
os.chdir(DATA)

dataset = sys.argv[1]
os.chdir(dataset)

OUTPUT = DATA + dataset + '/crystalfifps/fuzzy_ifp.fp'

processed = [o for o in os.listdir("./processed/") if os.path.isfile(os.path.join("./processed/",o))]
ligands = [o for o in os.listdir("./ligands/") if os.path.isfile(os.path.join("./ligands/",o))]

if not os.path.exists('crystalfifps'):
    os.system("mkdir crystalfifps")
os.chdir("crystalfifps")

for structFile in processed:
    struct = os.path.splitext(structFile)[0]
    ligandFile = filter(lambda x: struct.lower() in x.lower(), ligands)[0]
    #print(struct)
    #Change file below to: /scratch/PI/rondror/docking/workflow_thomas/fingerprint/fuzzy_ifp/fuzzyifp.py after I finish fixes
    with open(struct + '_fifpcrystaljob.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('OUTPUT=$(' + SCHRODINGER + ' ' + FUZZY_SCRIPT + ' -receptor ' + DATA + dataset + "/processed/" + structFile + ' -ligand ' + DATA + dataset + "/ligands/" + ligandFile + ' -input pair)\n')
        f.write("echo \"" + struct + ";" + '$OUTPUT\" >> ' + OUTPUT)
        
    os.system("sbatch --time=10:00:00 -n 1 -p rondror " + struct + '_fifpcrystaljob.sh')
    time.sleep(.5)
