import os
from multiprocessing import Pool
import subprocess

from schrodinger.structure import StructureReader
from schrodinger.application.prepwizard import fix_common_structure_mistakes

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(struct):
    subprocess.call([os.path.join(SCHRODINGER, "utilities", "prepwizard"),
                     "-WAIT", "-j", "temp-{}".format(struct),
                     "../aligned_proteins/{}.mae".format(struct),
                     "{}.mae".format(struct)])

    if '{}.mae'.format(struct) in os.listdir('.'):
        st = StructureReader('{}.mae'.format(struct)).next()
        if st._getTitle() != struct:
            os.system('rm {}.mae'.format(struct))

    return struct, '{}.mae'.format(struct) in os.listdir('.')

def process():
    os.system('mkdir -p processed_proteins')

    #for f in os.listdir('aligned_proteins'):
    #    if f in os.listdir('processed_proteins'):
    #        continue
    #    st = StructureReader('aligned_proteins/{}'.format(f)).next()
    #    print f, fix_common_structure_mistakes(st)
        #break

    os.chdir('processed_proteins')

    unfinished_structs = [f.split('.')[0] for f in os.listdir('../aligned_proteins') if f not in os.listdir('.')]
    pool = Pool(int(os.environ.get("SLURM_NTASKS", 4)))#, 10)

    for i in range(20):
        processing_structs = unfinished_structs
        unfinished_structs = []
        print('iteration {}, processing {} structs'.format(i+1, len(processing_structs)))
        print(processing_structs)
        for s, done in pool.imap_unordered(processHelper, processing_structs):
            print s, done
            if not done:
                unfinished_structs.append(s)
        if len(unfinished_structs) == 0:
            break
    else:
        print unfinished_structs, 'failed to process'

    os.system('rm temp* *missing* *.out')
    os.chdir('..')
