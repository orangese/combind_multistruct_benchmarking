#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
import time

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1/run"
FUZZY_SCRIPT = "/share/PI/rondror/$USER/combind/2_fingerprint/fuzzyifp.py"
DATA = "/scratch/PI/rondror/docking_data/"
if sys.argv[1] == 'all':
    datasets = os.listdir(DATA)
else:
    datasets = [sys.argv[1]]

output_dir = sys.argv[2]

for d in datasets:
    print d
    if d == 'B2AR_new': continue

    os.chdir(DATA + d)

    OUTPUT = DATA + d + '/' + output_dir + '/ifp.fp'

    processed = [o for o in os.listdir("./processed/") if os.path.isfile(os.path.join("./processed/",o))]
    ligands = [o for o in os.listdir("./ligands/") if os.path.isfile(os.path.join("./ligands/",o))]

    if os.path.exists(output_dir):
        #os.system('rm -r '+output_dir)
        print d + ' already has that output folder. continuing.'
        continue

    os.system("mkdir " + output_dir)
    os.chdir(output_dir)

    for structFile in processed:
        struct = os.path.splitext(structFile)[0]

        verbose_output = DATA + d + '/' + output_dir + '/' + struct + '.out'

        ligandFile = filter(lambda x: struct.lower() in x.lower(), ligands)[0]
        with open(struct + '.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('OUTPUT=$(' + SCHRODINGER + ' ' + FUZZY_SCRIPT + ' -receptor ' + DATA + d + "/processed/" + structFile + ' -ligand ' + DATA + d + "/ligands/" + ligandFile + ' -input pair -verbose ' + verbose_output + ')\n')
            f.write("echo \"" + struct + ";" + '$OUTPUT\" >> ' + OUTPUT)
        
        os.system("sbatch --time=00:10:00 -n 1 -p rondror " + struct + '.sh')


