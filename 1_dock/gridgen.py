#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os
import slurm
from multiprocessing import Pool

GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

reference_ligands = {
    'B1AR_all': '2Y00',
    'B2AR_all': '3PDS',
    'CHK1_all': '2C3K'
}

def getCentroid(receptor): 
    ref_ligand = '../ligands/{}_ligand.mae'.format(reference_ligands[receptor])
    struct = StructureReader(ref_ligand).next()
    asl_searcher = AslLigandSearcher()
    ligands = asl_searcher.search(struct)

    if len(ligands) == 0:
        print "Error: Could not find a reference ligand for {}".format(receptor)
        raise Exception()

    ligand = ligands[0].mol_num

    position_sum = [0, 0, 0]
    for atom in struct.molecule[ligand].atom:
        position_sum = [i+j for i, j in zip(position_sum, atom.xyz)]

    return map(lambda x: x / float(len(struct.molecule[ligand].atom)), position_sum)

def generateInFiles(structs, receptor):

    x, y, z = getCentroid(receptor)

    for s in structs:
        if '{}.in'.format(s) in os.listdir('.'): continue

        out = open("{}.in".format(s), 'w')
        out.write('GRID_CENTER {},{},{}\n'.format(x,y,z))
        out.write('GRIDFILE {}.zip\n'.format(s))

        struct = StructureReader('../processed/{}.mae'.format(s)).next()
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)
        if len(ligands) > 1:
            print "Error: multiple ligands in {}".format(s)
            raise Exception()
        if len(ligands) == 1:
            out.write('LIGAND_MOLECULE {}\n'.format(ligands[0].mol_num))

        out.write('INNERBOX 15,15,15\n')
        out.write('OUTERBOX 40,40,40\n')        
        out.write('RECEP_FILE ../processed/{}.mae\n'.format(s))
        out.close()

def getGridsHelper(struct):
    if struct in os.listdir('.'):
        return struct, True

    slurm.salloc("{} -WAIT {}.in".format(GLIDE, struct), "1", "1:00:00")

    if '{}.zip'.format(struct) in os.listdir('.'):
        os.system('mkdir {}'.format(struct))
        os.system('mv {}.zip {}'.format(struct, struct))
        return struct, True
    return struct, False

def getGrids(receptor):
    os.system('mkdir -p grids')
    os.chdir('grids')

    unfinished_grids = [f.split('.')[0] for f in os.listdir('../processed') if f.split('.')[1] == 'mae']
    
    generateInFiles(unfinished_grids, receptor)
    pool = Pool(min(len(unfinished_grids), 25))
    
    for i in range(1):
        processing_grids = unfinished_grids
        unfinished_grids = []
        print 'iteration {}, generating {} grids'.format(i+1, len(processing_grids))
        print processing_grids

        for g, done in pool.imap_unordered(getGridsHelper, processing_grids):
            if done: print g, 'succeeded!'
            else: unfinished_grids.append(g)
        if len(unfinished_grids) == 0:
            print 'all done!'
            os.system('rm *.log *.in gpu*')
            break
    else:
        print unfinished_grids, 'failed to generate grids'
    
    os.chdir('..')
