#!/share/PI/rondror/software/schrodinger2017-1/run

import sys
import os

PDBBIND = '/scratch/PI/rondror/docking_data/PDBbind_refined' # 2016 release
uniprot_index = '{}/index/INDEX_refined_name.2016'.format(PDBBIND)

command = sys.argv[1]

output = '/scratch/PI/rondror/docking_data/pdbbind'

os.chdir(PDBBIND)

def parse_index_file():
    test = {}
    with open(uniprot_index) as f:
        for line in f:
            if line[0] == '#': continue
            pdb_code, year, uniprot_id = line.split()[:3]
            if uniprot_id == '------': 
                uniprot_id = 'misc'
            if uniprot_id not in test:
                test[uniprot_id] = []
            test[uniprot_id].append(pdb_code)
    return test

os.chdir(output)

if command == '-m':
    for s in os.listdir('.'):
        os.system('rm -r {}/ligands {}/raw_pdbs'.format(s, s))
        os.system('mkdir {}/ligands'.format(s))
        os.system('mkdir {}/raw_pdbs'.format(s))
        for pdb in os.listdir('{}/original_files'.format(s)):
            os.system('cp {}/original_files/{}/{}_ligand.mol2 {}/ligands/{}_ligand.mol2'.format(s, pdb, pdb, s, pdb.upper()))
            os.system('cp {}/original_files/{}/{}_protein.pdb {}/raw_pdbs/{}.pdb'.format(s, pdb, pdb, s, pdb.upper()))
if command in ['-r', '-s']:
    num_folders = len(os.listdir('.'))
    for i, s in enumerate(os.listdir('.')):
        num_structs = len(os.listdir('{}/raw_pdbs'.format(s)))
        print 'running {} on {} structures in dataset {}. {} datasets remaining.'.format(command, num_structs, s, num_folders - i)
        os.system('./main.py {} pdbbind/{}'.format(command, s))
