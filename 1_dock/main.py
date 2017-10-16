##!./schrodinger_wrapper.sh

import sys
import os

from multiprocessing import Pool

import get_pdbs
import align_structures
import ligprep
import gridgen
import processing
import renumber
import dock
import tests
import prep_pdbbind
import stats
import protein_rmsds
import unique_ligands
import fp_2d

SCHRODINGER = os.environ.get("SCHRODINGER", None)
DATA = "/scratch/PI/rondror/docking_data"

command = sys.argv[1]
datasets = sys.argv[2:]

if datasets[0] == 'pdbbind_combo':
    datasets = ['pdbbind_combo/{}'.format(s) for s in os.listdir('{}/pdbbind_combo'.format(DATA))]
    datasets.sort()

if datasets[0] == 'pdbbind_deck':
    #datasets = [s for s in os.listdir('{}/pdbbind'.format(DATA)) 
    #            if len(os.listdir('{}/pdbbind/{}/ligands'.format(DATA, s))) >= 15 
    #            and s not in os.listdir('{}/pdbbind_final'.format(DATA))]
    #os.system('mkdir -p {}/pdbbind_deck'.format(DATA))
    #for d in datasets:
    #    print d
    #    if d not in os.listdir('{}/pdbbind_deck'.format(DATA)):
    #        os.system('cp -r {}/pdbbind/{} {}/pdbbind_deck'.format(DATA, d, DATA))
    #    for i in os.listdir('{}/pdbbind_deck/{}'.format(DATA, d)):
    #        if i not in ['ligands','original_files','raw_maes','raw_pdbs','processed_ligands','broken_ligands']:
    #            print 'deleting', i
    #            os.system('rm -rf {}/pdbbind_deck/{}/{}'.format(DATA, d, i))
    #        #print i

    datasets = ['pdbbind_deck/{}'.format(s) for s in os.listdir('{}/pdbbind_deck'.format(DATA))]
    datasets.sort()
#for d in datasets:
#    print d
#    os.system('rm -rf {}/{}/final_proteins {}/{}/final_ligands'.format(DATA, d, DATA, d))
#    os.system('mkdir {}/{}/final_proteins {}/{}/final_ligands'.format(DATA, d, DATA, d))
#    os.system('cp {}/{}/processed_proteins/* {}/{}/final_proteins'.format(DATA, d, DATA, d))
#    os.system('cp {}/{}/unique_ligands/* {}/{}/final_ligands'.format(DATA, d, DATA, d))
#    for l in os.listdir('{}/{}/final_ligands'.format(DATA, d)):
#        if l not in os.listdir('{}/{}/unique_ligands'.format(DATA, d)):
#            print l
#            #os.system('rm {}/{}/final_ligands/{}'.format(DATA, d, l))

command_map = {
    'b' : prep_pdbbind.move_files,
    'ge': get_pdbs.get,
    'l' : ligprep.extract_ligands,
    'u' : unique_ligands.manage_ligands,
    'a' : align_structures.align,
    's' : stats.get_alignment_stats,
    'p' : processing.process,
    'cp': unique_ligands.process_proteins,
    'n' : renumber.renumber,
    'r' : protein_rmsds.calc_rmsds,
    'gf': tests.get_first,
    'g' : gridgen.get_grids,
    'x' : dock.dock_dataset,
    '2d': fp_2d.get_2d_fp
}

def run_command(dataset):
    print command, dataset
    os.chdir('{}/{}'.format(DATA, dataset))    
    #print len(os.listdir('unique_ligands'))
    #return
    command_map[command]()
    return dataset

if len(datasets) == 1:
    run_command(datasets[0])
else:
    for d in datasets:
        if command == 'nvm': continue
        print d
        if command in ['x', 'u', 'cp', 'gf']:
            run_command(d)
        else:
            name = d
            if '/' in d:
                name = d.split('/')[1]
            os.system('sbatch --job-name={}-{} xdock.sh {} {}'.format(command, name, command, d))
