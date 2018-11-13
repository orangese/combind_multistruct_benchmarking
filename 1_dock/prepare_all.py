import os
import sys

from shared_paths import shared_paths, proteins
from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files

from grids import make_grids
from dock import dock, verify_dock
from fp_controller import compute_fp

from chembl_sort import get_ligands, proc_ligands
from chembl_props import write_props
from pick_helpers import pick_helpers, load_helpers

from containers import Protein

os.chdir(shared_paths['data'])

todo = sys.argv[1]

datasets = sys.argv[2:]
if datasets == []:
    datasets = proteins

for i, d in enumerate(datasets):
    print(d, i)
    os.chdir(d)
    protein = Protein(d, shared_paths['pdb_order'])
    lm = protein.lm

    if '0' in todo:
        sort_downloads()
 
    # 1. prepare proteins and associated grids 
    if '1' in todo:
        process_structs()             # Runs prepwizard
        align_structs()               # Align and give consistent numbering
        sort_files()                  # Creates ligand, protein, and complex directories
        make_grids()                  # Creates grid for all proteins
     
    # 2. prepare ligands
    if '2' in todo:
       get_ligands()                  # Writes MAE files for all ligs to ligands/raw_files
       proc_ligands()                 # Runs prepwizard & epik on all ligs
       
    if 'm' in todo:
        lm.mcss.compute_mcss() # Computes MCSS, for use in pick_helpers

    if 'p' in todo:
        dock(lm)
        dock(lm, mode = 'confgen_es1')
        dock(lm, mode = 'confgen_es4')
        dock(lm, mode = 'inplace')
        dock(lm, mode = 'mininplace')
        dock(lm, mode = 'XP')
        dock(lm, mode = 'expanded')
        lm.mcss.compute_mcss(False)
        compute_fp(lm, raw = 'raw' in shared_paths['ifp']['version'])

    if 'v' in todo:
        verify_mcss(lm)
    if 'g' in todo:
        check_docked_ligands(lm)
    # force redo of chembl info (do this if new chembl ligands have been added)
    # Do this after all MCSS files have been written!
    if 'c' in todo: #pass
         os.system('rm chembl/helpers/*')
         os.system('rm chembl/duplicates.txt')
         os.system('rm chembl/molw.txt')
         os.system('rm chembl/macrocycle.txt') 
         write_props(lm)

    # # 3. decide what ligands to use and prepare them
    if '3' in todo:
        pick_helpers(lm)                  # Picks chembl ligands for use in scoring for each pdb ligand
        dock(lm, load_helpers())          # Dock chembl ligands to be used for scoring all pdb ligands
        compute_fp(lm)                    # Writeout fingerprints for docked poses and pdb structures
        lm.mcss.compute_mcss(load_helpers())  # Performs all phases of MCSS computation

    os.chdir('..')
