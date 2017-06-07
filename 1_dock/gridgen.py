#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os

GLIDE = "/share/PI/rondror/software/schrodinger2016-1/glide"

def generateIn(struct_bases, pool):
    pool.map(generateInHelper, struct_bases)

def generateInHelper(struct_base):
    struct = StructureReader(struct_base+'.mae').next()

    '''
    # Identify signal ligand
    ligand = None
    for mol in struct.molecule:
        if 5 < len(filter(lambda x: x.element != 'H', mol.atom)) < 100:
            assert not ligand, 'There are two ligand sized molecules!'
            ligand = mol.number
    if not ligand:
        return False
    '''

    asl_searcher = AslLigandSearcher()
    ligands = asl_searcher.search(struct)
    if(len(ligands) > 1):
        assert True, 'There are %r ligand sized molecules' % len(ligands)
    if(len(ligands) == 0):
        assert True, 'Could not find a ligand'

    ligand = ligands[0].mol_num

    position_sum = [0, 0, 0]
    for atom in struct.molecule[ligand].atom:
        position_sum = [i+j for i, j in zip(position_sum, atom.xyz)]
    center = map(lambda x: x / float(len(struct.molecule[ligand].atom)), position_sum)

    out = open("{}.in".format(struct_base), 'w')
    out.write('GRID_CENTER '     + ','.join(map(str, center)) + '\n')
    out.write('GRIDFILE '        + struct_base                + '.zip\n')
    out.write('LIGAND_MOLECULE ' + str(ligand)                + '\n')
    out.write('RECEP_FILE '      + struct_base                + '.mae\n')
    out.close()

    return True

def runGlide(inFile):
    os.system(GLIDE + " -WAIT " + inFile)
