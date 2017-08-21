import os
from multiprocessing import Pool
import subprocess

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(struct):
    subprocess.call([os.path.join(SCHRODINGER, "utilities", "prepwizard"),
                     "-WAIT", "-fillsidechains", "-f", "3", "-fix",
                     "-samplewater", "-delwater_hbond_cutoff", "2",
                     "-keepfarwat", "-captermini",
                     "-j", "temp-{}".format(struct),
                     "../aligned_proteins/{}.mae".format(struct),
                     "{}.mae".format(struct)
                     ])

    return struct, '{}.mae'.format(struct) in os.listdir('.')

def process():
    os.system('mkdir -p processed')
    os.chdir('processed')

    unfinished_structs = [f.split('.')[0] for f in os.listdir('../aligned_proteins') if f not in os.listdir('.')]
    pool = Pool(int(os.environ.get("SLURM_NTASKS", 4)))#, 10)

    for i in range(2):
        processing_structs = unfinished_structs
        unfinished_structs = []
        print('iteration {}, processing {} structs'.format(i+1, len(processing_structs)))
        print(processing_structs)
        for s, done in pool.imap_unordered(processHelper, processing_structs):
            if not done:
                unfinished_structs.append(s)
        if len(unfinished_structs) == 0:
            break
    else:
        print unfinished_structs, 'failed to process'

    os.system('rm temp* *missing* *.out')
    os.chdir('..')
