#!./schrodinger_wrapper.sh

import sys
import os

from multiprocessing import Pool

import get_pdbs
import align_structures
import ligprep
import gridgen
import processing
import dock
import tests
import prep_pdbbind
import stats

SCHRODINGER = os.environ.get("SCHRODINGER", None)
DATA = "/scratch/PI/rondror/docking_data"

command = sys.argv[1][1]
datasets = sys.argv[2:]

if datasets[0] == 'pdbbind_final':
    datasets = ['pdbbind_final/{}'.format(s) for s in os.listdir('{}/pdbbind_final'.format(DATA))]

print datasets

command_map = {
    'm' : prep_pdbbind.move_files,
    'r' : get_pdbs.get,
    'l' : ligprep.extract_ligands,
    'a' : align_structures.align,
    's' : stats.get_alignment_stats,
    'p' : processing.process,
    'g' : gridgen.get_grids,
    'x' : dock.dock_dataset
}

all_commands = ['m','r','l','a','p','s','g','x']

def run_command(dataset):
    os.chdir('{}/{}'.format(DATA, dataset))    
    command_map[command]()
    return dataset

pool = Pool(min(10,len(datasets)))
i=0
for d in pool.imap_unordered(run_command, datasets):
    i+=1
    print i, d

