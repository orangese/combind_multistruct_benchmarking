#!/share/PI/rondror/software/schrodinger2017-1/run

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

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
DATA = "/scratch/PI/rondror/docking_data"

command = sys.argv[1][1]
datasets = sys.argv[2:]

if datasets[0] == 'pdbbind':
    all_pdbbind = tests.parse_index_file()
    datasets = ['pdbbind/{}'.format(s) for s in os.listdir('{}/pdbbind'.format(DATA))]
    datasets = [d for d in datasets if 'aligned_ligands' in os.listdir('{}/{}'.format(DATA,d))]
    datasets = [d for d in datasets if len(os.listdir('{}/{}/aligned_ligands'.format(DATA,d))) >= 8]
    print len(datasets)

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
    
    already_done = tests.check_prerequisites(all_commands[all_commands.index(command) + 1])
    if already_done: return
    if not tests.check_prerequisites(command):
        print '{} not ready to {}'.format(dataset, command)
        return
    print 'running {} on {}'.format(command, dataset)
    try:
        command_map[command]()
    except Exception as e:
        print e

for i, d in enumerate(datasets):
    if d == 'pdbbind/misc': continue
    print d, len(datasets) - i
    run_command(d)
#pool = Pool(1)#int(os.environ.get("SLURM_NTASKS",4))) # , 10)
#pool.imap_unordered(run_command, datasets)
