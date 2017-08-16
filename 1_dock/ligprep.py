import os
import sys
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter

# New goals:
# 1. Make sure proteins are aligned
# 2. Make sure all ligands are in the orthosteric site

# Plan:
# compare all detected ligands (1 or multiple) to the reference ligand defined elsewhere
# ref. ligand is in gridgen.py right now. move it either to a different class, imported everywhere, or pass it from main. 

def extractLigands():

    os.system("mkdir -p ligands")
    structs = [f.split('.')[0] for f in os.listdir("processed") if not os.path.exists('ligands/{}_ligand.mae'.format(f.split('.')[0]))]

    print 'Extracting ligands from {} files'.format(len(structs))

    for s in structs:
        struct = StructureReader("processed/{}.mae".format(s)).next()
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)

        if len(ligands) == 0: print 'Error: could not find a ligand for {}'.format(s)

        elif len(ligands) > 1: print 'Error: multiple ligands found for {}'.format(s)
        
        else:
            print 'Found 1 ligand for {}. Success!'.format(s)
            ligand = ligands[0]
            st_writer = StructureWriter("ligands/{}_ligand.mae".format(s))
            st_writer.append(ligand.st)
            st_writer.close()
