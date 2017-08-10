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
    if len(ligands) == 0:
        print "Error: Could not find a ligand in {}".format(struct_base)
        return

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

def processSuccess(structure):
    return structure+".zip" in os.listdir(".") #We're currently in grids

def runGlidesHelper(inFile):
    slurm.salloc("{} -WAIT {}".format(GLIDE, inFile), "1", "1:00:00")
    return inFile.split(".")[0] #Return the sturcture name (before the .in extension)

def runGlides(inFiles):
    #Each thread has a very lightweight job (waiting on the salloc signal) so just launch len(inFiles) threads
    pool = Pool(len(inFiles))
    
    #Keep submitting jobs until they're all successful
    while len(inFiles) != 0:
        currentlyProcessing = inFiles
        inFiles = []

        for finishedStruct in pool.imap_unordered(runGlidesHelper, currentlyProcessing):
            if not processSuccess(finishedStruct):
                inFiles.append(finishedStruct)
                print("{} did not gridgen successfully - resubmitting".format(finishedStruct))

