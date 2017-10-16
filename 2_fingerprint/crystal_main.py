#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
import time

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1/run"
FUZZY_SCRIPT = "/share/PI/rondror/$USER/combind/2_fingerprint/fuzzyifp.py"
DATA = "/scratch/PI/rondror/docking_data/"

output_dir = sys.argv[1]
datasets = sys.argv[2:]

if datasets[0] == 'pdbbind_combo':
    datasets = ['pdbbind_combo/{}'.format(s) for s in os.listdir('{}/pdbbind_combo'.format(DATA))]
if datasets[0] == 'pdbbind_deck':
    datasets = ['pdbbind_deck/{}'.format(s) for s in os.listdir('{}/pdbbind_deck'.format(DATA))]
datasets.sort()


for d in datasets:
    print d

    os.system('mkdir -p {}{}/ifp'.format(DATA, d))
    os.chdir('{}{}/ifp'.format(DATA, d))
    if output_dir in os.listdir('.'):
        print d + ' already has that output folder. continuing.'
        #continue
        os.system('rm -r '+output_dir)
    os.system('mkdir ' + output_dir)
    os.chdir(output_dir)

    #assert len(os.listdir('../../grids')) == 1
    #grid = os.listdir('../../grids')[0]
    for lig in os.listdir('../../final_ligands'):
    #for structFile in os.listdir('../../processed'): 
    #    struct = os.path.splitext(structFile)[0]
    #    ligs = filter(lambda x: struct in x, os.listdir('../../aligned_ligands'))
    #    if ligs == []: continue
    #    lig = ligs[0]
        #if os.path.exists('{}{}/renumbered_proteins'.format(DATA, d)):
        struct_path = '{}{}/final_proteins/{}.mae'.format(DATA, d, lig.split('_')[0])
        #else:
        #    struct_path = '{}{}/processed_proteins/{}.mae'.format(DATA, d, lig.split('_')[0])

        lig_path = '{}{}/final_ligands/{}'.format(DATA, d, lig)

        with open('{}.sh'.format(lig.split('_')[0]), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('OUTPUT=$({} {} -receptor {} -ligand {} -input pair -verbose {}.out)\n'.format(SCHRODINGER, FUZZY_SCRIPT, struct_path, lig_path, lig.split('_')[0])) 
            f.write('echo \"$OUTPUT\" >> {}.fp'.format(lig.split('_')[0]))
        os.system('sbatch --time=01:00:00 -n 1 -p rondror {}.sh'.format(lig.split('_')[0]))
        #break

