import os
import sys

from dock.align_structs import align_structs
from dock.process_structs import process_structs
from dock.sort_files import sort_files

from dock.ligands import prep_ligands

from dock.grids import make_grids
from dock.dock import dock

from settings import paths, stats, proteins
from ifp.fp_controller import compute_fp
from containers import Protein

def main(args):
    stats_version = args[0]
    task = args[1]
    datasets = args[2:]
    if datasets == []:
        datasets = proteins
    params = stats[stats_version]
    
    for i, d in enumerate(datasets):
        print(d, i)
        os.chdir('{}/{}'.format(paths['DATA'], d))
        protein = Protein(d, params, paths)
     
        if task == 'prep-structs':
            process_structs()
            align_structs()
            sort_files()
            make_grids()

        if task == 'prep-ligands':
            prep_ligands(protein.lm)

        if task == 'mcss':
            protein.lm.mcss.compute_mcss()

        if task == 'pick-helpers':
            protein.lm.pick_helpers()

        if task == 'dock-pdb':
            dock(protein.lm, mode=params['docking_version'])
            protein.lm.mcss.compute_mcss(False)
            compute_fp(protein.lm)

        if task == 'dock-helpers':
            dock(protein.lm, protein.lm.load_helpers(), mode=params['docking_version'])
            protein.lm.mcss.compute_mcss(True, protein.lm.load_helpers())
            compute_fp(protein.lm)
