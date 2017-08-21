#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
import time

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1/run"
FUZZY_SCRIPT = "/share/PI/rondror/$USER/combind/2_fingerprint/fuzzyifp.py"
DATA = "/scratch/PI/rondror/docking_data/"

output_dir = sys.argv[1]
datasets = sys.argv[2:]

for d in datasets:
    print d

    os.system('mkdir -p {}{}/ifp'.format(DATA, d))
    os.chdir('{}{}/ifp'.format(DATA, d))
    if output_dir in os.listdir('.'):
        print d + ' already has that output folder. continuing.'
        continue

    os.system('mkdir ' + output_dir)
    os.chdir(output_dir)

    for structFile in os.listdir('../../processed'): 
        struct = os.path.splitext(structFile)[0]
        ligs = filter(lambda x: struct in x, os.listdir('../../aligned_ligands'))
        if ligs == []: continue
        lig = ligs[0]
        
        struct_path = '{}{}/processed/{}.mae'.format(DATA, d, struct)
        lig_path = '{}{}/aligned_ligands/{}'.format(DATA, d, lig)

        with open('{}.sh'.format(struct), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('OUTPUT=$({} {} -receptor {} -ligand {} -input pair -verbose {}.out)\n'.format(SCHRODINGER, FUZZY_SCRIPT, struct_path, lig_path, struct)) 
            f.write('echo \"{};$OUTPUT\" >> {}.fp'.format(struct, struct))
        os.system('sbatch --time=00:10:00 -n 1 -p rondror {}.sh'.format(struct))


