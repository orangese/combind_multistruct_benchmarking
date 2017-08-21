import os
import sys
import numpy as np
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter

def extract_ligands():

    os.system("mkdir -p ligands")
    structs = [f.split('.')[0] for f in os.listdir("raw_maes") if not os.path.exists('ligands/{}_ligand.mae'.format(f.split('.')[0]))]

    print 'Extracting ligands from {} files'.format(len(structs))
    for s in structs:
        struct = StructureReader("raw_maes/{}.mae".format(s)).next()
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)

        if len(ligands) == 0: 
	   print 'Error: could not find a ligand for {}'.format(s)
	   out = open("multiple_ligands.out", "a+")
	   out.write("Structure: {}".format(s))
	   out.write(" Number of ligands: 0")
	   out.write("\n")
	   out.close()

	elif len(ligands) == 1:
	   struct.deleteAtoms(ligands[0].atom_indexes)
	   struct.write("raw_maes/{}.mae".format(s))
	   st_writer = StructureWriter("ligands/{}_ligand.mae".format(s))
	   st_writer.append(ligands[0].st)
	   st_writer.close()

        else:
	     print "Multiple ligands detected in {}. See output file and manually remove them.".format(s)
	     out = open("multiple_ligands.out", "a+")
	     out.write("Structure: {}".format(s))
	     out.write(" Number of ligands: {}".format(len(ligands))) 
	     for ligand in ligands:
		out.write(" {}".format(ligand.pdbres))
	     out.write("\n")
	     out.close()
