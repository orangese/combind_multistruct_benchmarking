from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
import os
import subprocess
from multiprocessing import Pool
from tests import get_first
from centroid import average_atom_pos
GLIDE = '{}/glide'.format(os.environ.get("SCHRODINGER", None))# "/share/PI/rondror/software/schrodinger2017-1/glide"

def generate_input_files(structs, centroid):

    x, y, z = centroid

    for s in structs:
        if '{}.in'.format(s) in os.listdir('.'): continue

        out = open("{}.in".format(s), 'w')
        out.write('GRID_CENTER {},{},{}\n'.format(x,y,z))
        out.write('GRIDFILE {}.zip\n'.format(s))
        out.write('INNERBOX 15,15,15\n')
        out.write('OUTERBOX 40,40,40\n')
        out.write('RECEP_FILE ../final_proteins/{}.mae\n'.format(s))
        out.close()

def get_grids_helper(struct):
    if struct in os.listdir('.'):
        return struct, True

    st = StructureReader('../final_ligands/{}_ligand.mae'.format(struct)).next()
    centroid = average_atom_pos(st) 
    print struct, centroid
    generate_input_files([struct], centroid)

    subprocess.call([GLIDE, "-WAIT", "{}.in".format(struct)])

    if '{}.zip'.format(struct) in os.listdir('.'):
        os.system('mkdir {}'.format(struct))
        os.system('mv {}.zip {}'.format(struct, struct))
        return struct, True
    return struct, False

def get_grids(print_only=False):
    unfinished = [i.split('.')[0] for i in os.listdir('final_proteins') if i.split('.')[0] not in os.listdir('grids')]
    if os.path.exists('../../pdbbind_first_grids.txt'):
        grids1 = {}
        with open('../../pdbbind_first_grids.txt') as f:
            for line in f:
                dataset, grid = line.strip().split()
                grids1[dataset] = grid

        grids2 = {}
        with open('../../pdbbind_second_grids.txt') as f:
            for line in f:
                dataset, grid = line.strip().split()
                grids2[dataset] = grid

        dataset_name = os.getcwd().split('/')[-1]
        unfinished = [grids1[dataset_name],grids2[dataset_name]] 
    
    print unfinished
    if print_only: return
    os.system('mkdir -p grids')
    os.chdir('grids')

    pool = Pool(int(os.environ.get("SLURM_NTASKS", 2)))

    for i in range(2):
        processing_grids = unfinished
        unfinished = []
        print 'iteration {}, generating {} grids'.format(i+1, len(processing_grids))
        for g, done in pool.imap_unordered(get_grids_helper, processing_grids):
            if not done: 
                print g, 'failed'
                unfinished.append(g)
            else: 
                print g, 'succeeded'
        if len(unfinished) == 0:
            os.system('rm *.log *.in gpu* sherlock*')
            break
    else:
        print unfinished, 'failed to generate grids'

    os.chdir('..')
