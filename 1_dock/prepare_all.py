import os
import sys

from shared_paths import shared_paths
from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files

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

todo = list(sys.argv[1])
if len(todo) == 0:
    todo = list('12345') 

datasets = sys.argv[2:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir('.')) if d[0] != '.' and d[-3:] != 'old']

datasets=reversed(datasets)
for i, d in enumerate(datasets):
    print d, i
    os.chdir(d)
    
    lm = LigandManager(shared_paths, d)
 
    # 1. prepare proteins and associated grids 
    if '1' in todo: 
        pdb_st = None
        if os.path.exists('structures/downloads'):    
            pdb_st = sort_downloads()
        process_structs()             # Runs prepwizard
        align_structs()               # Align and give consistent numbering
        sort_files()                  # Creates ligand, protein, and complex directories
        make_grids()                  # Creates grid for all proteins
        dock(lm)                      # Dock PDB ligands to lm.st
     
    # 2. prepare ligands
    if '2' in todo:
       get_ligands()                  # Writes MAE files for all ligs to ligands/raw_files
       proc_ligands()                 # Runs prepwizard & epik on all ligs

    # force redo of chembl info (do this if new chembl ligands have been added)
    if 'c' in todo: #pass
         os.system('rm -f chembl/helpers/*')
         os.system('rm -f chembl/duplicates.txt')
         os.system('rm -f chembl/molw.txt')
         write_props(lm)

    # # 3. decide what ligands to use and prepare thems
    if '3' in todo:
        pick_helpers(lm)              # Picks chembl ligands for use in scoring for each pdb ligand
        dock(lm, load_helpers())      # Dock chembl ligands to be used for scoring all pdb ligands
        fp(lm)                        # Writeout fingerprints for docked poses and pdb structures
        mcss(lm, load_helpers())      # Performs all 3 phases of MCSS computation

    # 4. compute statistics
    if '4' in todo:
        stats(lm)

    # 5. run method on all ligands
    if '5' in todo:
        score(lm, load_helpers())

    os.chdir('..')
