#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os
import slurm
from multiprocessing import Pool

GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

def generateIn(struct_bases):
    #Single threaded since this doesn't take very long
    map(generateInHelper, struct_bases)

def generateInHelper(struct_base):
    struct = StructureReader(struct_base+'.mae').next()
    asl_searcher = AslLigandSearcher()
    ligands = asl_searcher.search(struct)
    assert len(ligands) != 0, "Error: Could not find a ligand in {}".format(struct_base)
    if(len(ligands) > 1):
        print('Warning: There are {} ligand sized molecules in {}, picking the first'.format(struct_base, len(ligands)))

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

def runGlidesHelper(inFile):
    slurm.salloc("{} -WAIT {}".format(GLIDE, inFile), "1", "1:00:00")

def runGlides(inFiles):
    #Each thread has a very lightweight job (waiting on the salloc signal) so just launch len(inFiles) threads
    pool = Pool(len(inFiles))
    pool.map(runGlidesHelper, inFiles)
