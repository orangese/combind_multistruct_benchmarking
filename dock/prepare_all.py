import os
import sys

from dock.sort_downloads import sort_downloads

from dock.align_structs import align_structs
from dock.process_structs import process_structs
from dock.sort_files import sort_files

from dock.grids import make_grids
from dock.dock import dock, verify_dock

from dock.chembl_sort import get_ligands, proc_ligands
from dock.chembl_props import write_props
from dock.pick_helpers import pick_helpers, load_helpers

from shared_paths import shared_paths, proteins
from ifp.fp_controller import compute_fp
from containers import Protein

def main(args):
    task = args[1]
    datasets = args[2:]
    if datasets == []:
        datasets = proteins
    os.chdir(shared_paths['data'])

    for i, d in enumerate(datasets):
        print(d, i)
        os.chdir(d)
        protein = Protein(d, shared_paths['pdb_order'])
        lm = protein.lm

        if task == '0':
            sort_downloads()
     
        if task == '1':
            process_structs()      # Runs prepwizard
            align_structs()        # Align and give consistent numbering
            sort_files()           # Creates ligand, protein, and complex directories
            make_grids()           # Creates grid for all proteins
         
        if task == '2':
            get_ligands()           # Writes MAE files for all ligs to ligands/raw_files
            proc_ligands()          # Runs prepwizard & epik on all ligs
           
        if task == 'm':
            lm.mcss.compute_mcss() # Computes MCSS, for use in pick_helpers

        if task == 'p':
            #dock(lm)
            #dock(lm, mode = 'confgen_es1')
            dock(lm, mode = 'confgen_es4')
            #dock(lm, mode = 'confgen_es4_soft')
            #dock(lm, mode = 'confgen_es4_nocanon')
            #dock(lm, mode = 'inplace')
            #dock(lm, mode = 'mininplace')
            dock(lm, mode = 'XP')
            #dock(lm, mode = 'expanded')
            lm.mcss.compute_mcss(False)
            compute_fp(lm, raw = 'raw' in shared_paths['ifp']['version'])

        if task == 'v':
            verify_mcss(lm)
        
        if task == 'g':
            check_docked_ligands(lm)
        
        # force redo of chembl info (do this if new chembl ligands have been added)
        # Do this after all MCSS files have been written!
        if task == 'c':
             os.system('rm chembl/helpers/*')
             os.system('rm chembl/duplicates.txt')
             os.system('rm chembl/molw.txt')
             os.system('rm chembl/macrocycle.txt') 
             write_props(lm)

        # 3. decide what ligands to use and prepare them
        if task == '3':
            pick_helpers(lm)         # Picks chembl ligands for use in scoring
            dock(lm, load_helpers()) # Dock chembl ligands to be used in scoring
            compute_fp(lm)           # Fingerprints for docked ligands and pdb structures
            lm.mcss.compute_mcss(True, load_helpers())
        
        # 3m. same as 3, except we have multiple mutant receptors to dock to
        if task == '3m':
            pick_helpers(lm)
            dock(lm, mutants=True)
            compute_fp(lm)
            lm.mcss.compute_mcss(True)

        os.chdir('..')
