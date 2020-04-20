import os
import sys

from dock.struct_align import struct_align
from dock.struct_process import struct_process
from dock.struct_sort import struct_sort

from dock.ligands import prep_ligands

from dock.grids import make_grids
from dock.dock import dock

from ifp.ifp_controller import compute_ifp
from containers import Protein

def main(params, paths, task, proteins):
    for i, protein_name in enumerate(proteins):
        print(protein_name, i)
        os.chdir('{}/{}'.format(paths['DATA'], protein_name))
        protein = Protein(protein_name, params, paths)
     
        if task == 'prep-structs':
            struct_process()
            struct_align()
            struct_sort()
            make_grids()

        elif task == 'prep-ligands':
            prep_ligands(protein.lm)

        elif task == 'dock':
            dock(protein.lm, mode=params['docking_version'])
            protein.lm.mcss.compute_mcss(False)
            compute_ifp(protein.lm)

        elif task == 'mcss':
            protein.lm.mcss.compute_mcss()

        elif task == 'pick-helpers':
            protein.lm.pick_helpers()

        elif task == 'dock-helpers':
            dock(protein.lm, protein.lm.load_helpers(), mode=params['docking_version'])
            protein.lm.mcss.compute_mcss(True, protein.lm.load_helpers())
            compute_ifp(protein.lm)
        else:
            assert False
