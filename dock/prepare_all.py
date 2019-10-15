import os
import sys

from dock.sort_downloads import sort_downloads

from dock.align_structs import align_structs
from dock.process_structs import process_structs
from dock.sort_files import sort_files

from dock.grids import make_grids
from dock.dock import dock

import dock.chembl_sort as chembl_sort
from dock.chembl_props import write_props

from settings import paths, stats, proteins
from ifp.fp_controller import compute_fp
from containers import Protein

def main(args):
    stats_version = args[1]
    task = args[2]
    datasets = args[3:]
    if datasets == []:
        datasets = proteins
    params = stats[stats_version]
    
    for i, d in enumerate(datasets):
        print(d, i)
        os.chdir('{}/{}'.format(paths['data'], d))
        protein = Protein(d, params['pdb_order'])

        if task == '0':
            sort_downloads()
     
        if task == '1':
            process_structs()      # Runs prepwizard
            align_structs()        # Align and give consistent numbering
            sort_files()           # Creates ligand, protein, and complex directories
            make_grids()           # Creates grid for all proteins
         
        if task == '2':
            chembl_sort.get_ligands()           # Writes MAE files for all ligs to ligands/raw_files
            chembl_sort.proc_ligands()          # Runs prepwizard & epik on all ligs

        if task == '2chembl':  # prep only chembl ligands from smiles
            chembl_sort.prep_chembl_workflow(paths['data']+'/'+d)
        if task == '2chembl_done_check':  # check if chembl prep done, print results
            chembl_sort.check_chembl_prep_complete(paths['data']+'/'+d)

        if task == 'm':
            lm.mcss.compute_mcss() # Computes MCSS, for use in pick_helpers

        if task == 'p':
            dock(protein.lm, mode=params['docking_version'])
            protein.lm.mcss.compute_mcss(False)
            compute_fp(protein.lm)
        
        # force redo of chembl info (do this if new chembl ligands have been added)
        # Do this after all MCSS files have been written!
        if task == 'c':
             os.system('rm chembl/helpers/*')
             os.system('rm chembl/duplicates.txt')
             os.system('rm chembl/molw.txt')
             os.system('rm chembl/macrocycle.txt') 
             write_props(protein.lm)

        # 3. decide what ligands to use and prepare them
        if task == '3':
            lm.pick_helpers()
            dock(protein.lm, protein.lm.load_helpers(),
                 mode=params['docking_version'])
            compute_fp(protein.lm)
            lm.mcss.compute_mcss(True, protein.lm.load_helpers())
        
        # 3m. same as 3, except we have multiple mutant receptors to dock to
        if task == '3m':
            pick_helpers(lm)
            dock(lm, mutants=True)
            compute_fp(lm)
            lm.mcss.compute_mcss(True)
