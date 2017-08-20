#!/share/PI/rondror/software/schrodinger2016-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os
import subprocess
from multiprocessing import Pool

GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

def generate_input_files(structs, centroid):

    x, y, z = centroid

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

    with open('grid_center.txt','r') as f:
        centroid = [float(i) for i in f[0].split(',')]

    os.system('mkdir -p grids')
    os.chdir('grids')

    unfinished_grids = [f.split('.')[0] for f in os.listdir('../processed') if f.split('.')[1] == 'mae']

    generate_input_files(unfinished_grids, centroid)
    pool = Pool(int(os.environ.get("SLURM_NTASKS", 4)))#, 10)

    for i in range(1):
        processing_grids = unfinished_grids
        unfinished_grids = []
        print 'iteration {}, generating {} grids'.format(i+1, len(processing_grids))
        print processing_grids

        for g, done in pool.imap_unordered(get_grids_helper, processing_grids):
            if done: print g, 'succeeded!'
            else: unfinished_grids.append(g)
        if len(unfinished_grids) == 0:
            print 'all done!'
            os.system('rm *.log *.in gpu*')
            break
    else:
        print unfinished_grids, 'failed to generate grids'

    os.chdir('..')
