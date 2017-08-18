import os
import sys
import numpy as np
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter


# New goals:
# 1. Make sure proteins are aligned
# 2. Make sure all ligands are in the orthosteric site

# Plan:
# compare all detected ligands (1 or multiple) to the reference ligand defined elsewhere
# ref. ligand is in gridgen.py right now. move it either to a different class, imported everywhere, or pass it from main. 

def extractLigands(ref_coord):

    os.system("mkdir -p ligands")
    structs = [f.split('.')[0] for f in os.listdir("processed") if not os.path.exists('ligands/{}_ligand.mae'.format(f.split('.')[0]))]

    print 'Extracting ligands from {} files'.format(len(structs))
    for s in structs:
        struct = StructureReader("processed/{}.mae".format(s)).next()
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)

        if len(ligands) == 0: print 'Error: could not find a ligand for {}'.format(s)

        else:
             (true_lig, distance) = (None, 5.0)
	     out = open("lig_distances.out", "a+")
	     out.write("Structure: {}".format(s))
	     out.write(" Number of ligands: {}".format(len(ligands))) 
	     for ligand in ligands:
		dist = np.linalg.norm(ref_coord - ligand.centroid) 
		out.write(" {}".format(ligand.pdbres))
		out.write(" {}".format(dist))
		if dist < distance:
		   (true_lig, distance) = (ligand, dist)
	     out.write("\n")
	     out.close()

	     if true_lig is None:
	        print 'Alignment for {} incorrect or no ligand in binding pocket. Align structures manually.'.format(s)
		#break    
             else:
                st_writer = StructureWriter("ligands/{}_ligand.mae".format(s))
       	        st_writer.append(true_lig.st)
                st_writer.close()
