import os
from multiprocessing import Pool
import subprocess

from schrodinger.structure import StructureReader

SCHRODINGER = os.environ.get("SCHRODINGER", None)
def process_helper(struct):
    subprocess.call([os.path.join(SCHRODINGER, "utilities", "prepwizard"), "-WAIT", "-fix",
                     "../raw_maes/{}.mae".format(struct), "{}.mae".format(struct)])

    if is_processed(struct):
        return struct, True
    
    #os.system('rm {}*'.format(struct))
    return struct, False

def is_processed(struct):
    if '{}.mae'.format(struct) not in os.listdir('.'): 
        print struct, 'file not present'
        return False
    
    st = StructureReader('{}.mae'.format(struct)).next()
    if st._getTitle() != struct: 
        os.system('rm {}.mae'.format(struct))
        print struct, 'wrong title'
        return False
    return True

def process(print_only=False):
    os.system('mkdir -p processed_proteins')
    os.chdir('processed_proteins')
    unfinished_structs = [f.split('.')[0] for f in os.listdir('../raw_maes') if not is_processed(f.split('.')[0])]
    if print_only:
        print len(unfinished_structs)
        return
    
    pool = Pool(int(os.environ.get("SLURM_NTASKS", 4)))
    
    for i in range(5):
        processing_structs = unfinished_structs
        unfinished_structs = []
        print('iteration {}, processing {} structs'.format(i+1, len(processing_structs)))
        print(processing_structs)
        for s, done in pool.imap_unordered(process_helper, processing_structs):
        #for struct in processing_structs:
            #s, done = process_helper(struct)
            print s, done
            if not done:
                unfinished_structs.append(s)
        if len(unfinished_structs) == 0:
            break
    else:
        print unfinished_structs, 'failed to process'

    for i in os.listdir('.'):
        if i not in set(os.listdir('../raw_maes')):
            os.system('rm {}'.format(i))

    os.chdir('..')
