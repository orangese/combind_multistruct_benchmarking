import os
import sys
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
from multiprocessing import Pool
import slurm

def extractLigands(structures):
    for struct_base in structures:
        indir = "processed"
        outdir = "ligands"
        struct = StructureReader("{}/{}.mae".format(indir, struct_base)).next()

        # Identify signal ligand
        ligand = None
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)

        assert len(ligands) != 0, 'Error: Could not find a ligand for {}'.format(struct_base)

        if(len(ligands) > 1):
            print('Warning: There are multiple ligand sized molecules, picking the first one')
        
        ligand = ligands[0]
        st_writer = StructureWriter("{}/{}_ligand.mae".format(outdir, struct_base))
        st_writer.append(ligand.st)
        st_writer.close()

