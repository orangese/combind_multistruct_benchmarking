from schrodinger.structure import StructureReader
import os

def check_docked_ligands(lm):
    for ligand in lm.docked(lm.all_ligs):
        prepared = "{0:}/ligands/prepared_ligands/{1:}/{1:}.mae".format(lm.root, ligand)
        docked = "{0:}/docking/{1:}/{2:}-to-{3:}/{2:}-to-{3:}_pv.maegz".format(lm.root,
                                                                               lm.sp['docking'],
                                                                               ligand,
                                                                               lm.st)
        
        if not os.path.exists(prepared) and os.path.exists(docked):
            print('Ligand or docking files do not exist for:')
            print(prepared)
            print(docked)
            continue

        prepared = next(StructureReader(prepared))
        docked = list(StructureReader(docked))[1]
        if not prepared.isEquivalent(docked):
            print('{} is not the same in glide output and prepared ligands.'.format(ligand))
