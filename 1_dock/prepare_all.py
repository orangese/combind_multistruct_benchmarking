import os
import sys

from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files
#from init_mcss import init_mcss

from grids import make_grids
from dock import dock
from fp import fp
from mcss import mcss
from stats import stats

from chembl_sort import get_ligands, proc_ligands
from chembl_props import write_props
from pick_helpers import pick_helpers, load_helpers
from score import score

sys.path.append('../3_analyze')
from containers import LigandManager

os.chdir('../../data')

shared_paths = {
    'code'    : '/scratch/PI/rondror/jbelk/method/combind',
    'data'    : '/scratch/PI/rondror/jbelk/method/data',
    'docking' : 'glide12',
    'mcss'    : 'custom_mcss',
    'stats'   : 'stats3',
    'ifp'     : 'ifp3'
}

todo = list(sys.argv[1])
if len(todo) == 0:
    todo = list('12345') 

datasets = sys.argv[2:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir('.')) if d[0] != '.' and d[-3:] != 'old']

for i, d in enumerate(datasets):
    print i, d
    os.chdir(d)
    
    lm = LigandManager(shared_paths, d)
 
    # 1. prepare proteins   
    if '1' in todo: 
        pdb_st = None
        if os.path.exists('structures/downloads'):    
            pdb_st = sort_downloads()
        process_structs()
        align_structs()
        sort_files()
        make_grids()
        dock(lm)
     
    # 2. prepare ligands
    if '2' in todo:
        get_ligands()    
        proc_ligands()
        #mcss(lm)

    # force redo of chembl info (do this if new chembl ligands have been added)
    if 'c' in todo: #pass
        os.system('rm -f chembl/helpers/*')
        #os.system('rm -f chembl/duplicates.txt')
        #os.system('rm -f chembl/molw.txt')
        #write_props(lm)

    # 3. decide what ligands to use and prepare them
    if '3' in todo:
        pick_helpers(lm)
        #init_mcss(lm, load_helpers())
        dock(lm, load_helpers())
        fp(lm)
        mcss(lm, load_helpers())    

    # 4. compute statistics
    if '4' in todo:
        stats(lm)

    # 5. run method on all ligands
    if '5' in todo:
        score(lm, load_helpers())

    os.chdir('..')











