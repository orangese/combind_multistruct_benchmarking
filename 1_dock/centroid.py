#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter

#Ideally the reference ligand should be the only ligand in the structure.
#Otherwise, we need to add code that specifies which ligand in the structure is the reference ligand.
reference_ligands = {
    'B1AR_all': '2Y00',
    'B2AR_all': '3PDS',
    'CHK1_all': '2C3K',
    'TRMD_all': '3AXZ',
    'TRMD': '3AXZ'
}

def getCentroid(receptor): 
    #ref_ligand = 'ligands_all/{}_ligand.mae'.format(reference_ligands[receptor])
    struct = StructureReader("processed/{}.mae".format(reference_ligands[receptor])).next()
    asl_searcher = AslLigandSearcher()
    ligands = asl_searcher.search(struct)

    if len(ligands) == 0:
        print "Error: Could not find a reference ligand for {}".format(receptor)
        raise Exception()

    ligand = ligands[0]
    return ligand.centroid
