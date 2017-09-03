#!/share/PI/rondror/software/schrodinger2017-1/run
##!/share/software/modules/chemistry/schrodinger/2017-2/run
import sys
import os

from multiprocessing import Pool

from schrodinger.structure import StructureReader, StructureWriter, SmilesReader, SmilesWriter
from schrodinger.structutils.measure import get_shortest_distance
import schrodinger.structutils.analyze as schro
import schrodinger.structutils.interactions.pi as pi

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

p1 = ['pdbbind_final/{}'.format(s) for s in os.listdir('{}/pdbbind_final'.format(DATA))]
p2 = ['pdbbind_part2/{}'.format(s) for s in os.listdir('{}/pdbbind_part2'.format(DATA))]
#print p1
os.chdir(DATA)
for i in p1 + p2:
    print i
    os.chdir(i)
    dock.dock_dataset(True)
    os.chdir(DATA)
    #pass 

#tests.check_duplicates('{}/pdbbind_final'.format(DATA))
#tests.check_duplicates('{}/pdbbind_part2'.format(DATA))

if datasets[0] == 'pdbbind_part2':
    datasets = ['pdbbind_part2/{}'.format(s) for s in reversed(os.listdir('{}/pdbbind_part2'.format(DATA)))]

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

for i, d in enumerate(datasets):
    print d, len(datasets) - i
    run_command(d)
