#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os
import subprocess
from multiprocessing import Pool
from tests import get_first
from centroid import average_atom_pos
GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

def generate_input_files(structs, centroid):

    x, y, z = centroid

    for s in structs:
        if '{}.in'.format(s) in os.listdir('.'): continue

        out = open("{}.in".format(s), 'w')
        out.write('GRID_CENTER {},{},{}\n'.format(x,y,z))
        out.write('GRIDFILE {}.zip\n'.format(s))

        struct = StructureReader('../processed_proteins/{}.mae'.format(s)).next()
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)
        if len(ligands) > 1:
            print "Error: multiple ligands in {}".format(s)
            raise Exception()
        if len(ligands) == 1:
            out.write('LIGAND_MOLECULE {}\n'.format(ligands[0].mol_num))

        out.write('INNERBOX 15,15,15\n')
        out.write('OUTERBOX 40,40,40\n')
        out.write('RECEP_FILE ../processed_proteins/{}.mae\n'.format(s))
        out.close()

def get_grids_helper(struct):
    if struct in os.listdir('.'):
        return struct, True

    subprocess.call([GLIDE, "-WAIT", "{}.in".format(struct)])

    if '{}.zip'.format(struct) in os.listdir('.'):
        os.system('mkdir {}'.format(struct))
        os.system('mv {}.zip {}'.format(struct, struct))
        return struct, True
    return struct, False

def get_grids():

    #with open('grid_center.txt','r') as f:
    #    for line in f:
    #        centroid = [float(i) for i in line.strip().split(',')]
    #        break

    #if all_grids:
    #    unfinished_grids = [f.split('.')[0] for f in os.listdir('processed_proteins') if f.endswith('mae')]
    #else:
    #g = get_first()
    #print g
    #os.system('cp aligned_proteins/{}.mae ../'.format(g))
    #os.system('rm -rf processed_proteins')
    #os.mkdir('processed_proteins')
    #os.system('cp ../{}.mae processed_proteins'.format(g))
    #return
    only_grid = os.listdir('processed_proteins')
    assert len(only_grid) == 1
    only_grid = only_grid[0].split('.')[0]
    #unfinished_grids = [get_first()]
    #assert unfinished_grids == [os.listdir('processed_proteins')[0].split('.')[0]]
    #print unfinished_grids
    st = StructureReader('processed_ligands/{}_ligand.mae'.format(only_grid)).next()
    centroid = average_atom_pos(st) 
    print centroid
    #if len(os.listdir('grids')) == 2:
    #    print os.listdir('grids')
    #    os.system('rm -rf grids/{}'.format(unfinished_grids[0]))
    #return
    os.system('mkdir -p grids')
    os.chdir('grids')
    generate_input_files([only_grid], centroid)
    print only_grid
    #pool = Pool(int(os.environ.get("SLURM_NTASKS", 4)))#, 10)

    #for i in range(2):
    #    processing_grids = unfinished_grids
    #    unfinished_grids = []
    #    print 'iteration {}, generating {} grids'.format(i+1, len(processing_grids))
    #    print processing_grids
    struct, done = get_grids_helper(only_grid)
    #    for g, done in pool.imap_unordered(get_grids_helper, processing_grids):
    if not done: 
        print 'failed'
    else: 
        print 'succeeeded'
        os.system('rm *.log *.in gpu* sherlock*')
    #            unfinished_grids.append(g)
    #    if len(unfinished_grids) == 0:
    #        os.system('rm *.log *.in gpu* sherlock*')
    #        break
    #else:
    #    print unfinished_grids, 'failed to generate grids'

    os.chdir('..')
