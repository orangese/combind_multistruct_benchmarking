import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter

def extract(struct_base):
    indir = "processed"
    outdir = "ligands"
    struct = StructureReader("{}/{}.mae".format(indir, struct_base)).next()

    # Identify signal ligand
    ligand = None
    '''
    for mol in struct.molecule:
        if 5 < len(filter(lambda x: x.element != 'H', mol.atom)) < 100:
            assert not ligand, 'There are two ligand sized molecules!'
            ligand = mol.extractStructure(copy_props=True)
    if not ligand:
        return False
    '''
    asl_searcher = AslLigandSearcher()
    ligands = asl_searcher.search(struct)
    if(len(ligands) > 1):
        assert True, 'There is more than 1 ligand sized molecules'
    if(len(ligands) == 0):
        assert True, 'Could not find a ligand'

    print(struct_base + " " + str(ligands))
    ligand = ligands[0]
    st_writer = StructureWriter("{}/{}_ligand.mae".format(outdir, struct_base))
    st_writer.append(ligand.st)
    st_writer.close()

    # ligand.write("{}/{}_ligand.mae".format(outdir, struct_base))

    return True
