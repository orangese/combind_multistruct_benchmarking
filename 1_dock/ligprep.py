import os
import sys
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter

def extractLigands():

    os.system("mkdir -p ligands")
    structs = [f.split('.')[0] for f in os.listdir("processed") if not os.path.exists('ligands/{}_ligand.mae'.format(f.split('.')[0]))]

    print('extracting ligands from: ',structs)

    for s in structs:
        struct = StructureReader("processed/{}.mae".format(s)).next()
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)

        assert len(ligands) != 0, 'Error: Could not find a ligand for {}'.format(s)

        if(len(ligands) > 1):
            print('Warning: There are multiple ligand sized molecules, picking the first one')
        
        ligand = ligands[0]
        st_writer = StructureWriter("ligands/{}_ligand.mae".format(outdir, s))
        st_writer.append(ligand.st)
        st_writer.close()

