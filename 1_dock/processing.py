import os
from multiprocessing import Pool
import subprocess

from schrodinger.structure import StructureReader

SCHRODINGER = os.environ.get("SCHRODINGER", None)
def process_helper(struct):
    subprocess.call([os.path.join(SCHRODINGER, "utilities", "prepwizard"), "-WAIT", "-fix",
                     "../aligned_proteins/{}.mae".format(struct), "{}.mae".format(struct)])

    if is_processed(struct):
        os.system('rm {}-protassign* {}-impref* {}.log'.format(struct, struct, struct))
        return struct, True
    
    os.system('rm {}*'.format(struct))
    return struct, False

def is_processed(struct):
    if '{}.mae'.format(struct) not in os.listdir('.'): return False
    
    st = StructureReader('{}.mae'.format(struct)).next()
    if st._getTitle() != struct: return False

    valence = {'H':set([1]), 'C':set([4]), 'N':set([3]), 'O':set([2]), 'S':set([2,4,6])}

    for a in st.atom:
        if a.element in valence and sum([b.order for b in a.bond]) - a.formal_charge not in valence[a.element]:
            print struct, a.element, a.bond_total, a.formal_charge, a.index
            return False
    return True 

def process(print_only=False):
    os.system('mkdir -p processed_proteins')
    os.chdir('processed_proteins')
    unfinished_structs = [f.split('.')[0] for f in os.listdir('../aligned_proteins') if not is_processed(f.split('.')[0])]
    if print_only:
        print len(unfinished_structs)
        return
    pool = Pool(10)
    
    for i in range(10):
        processing_structs = unfinished_structs
        unfinished_structs = []
        print('iteration {}, processing {} structs'.format(i+1, len(processing_structs)))
        print(processing_structs)
        for s, done in pool.imap_unordered(process_helper, processing_structs):
            print s, done
            if not done:
                unfinished_structs.append(s)
        if len(unfinished_structs) == 0:
            break
    else:
        print unfinished_structs, 'failed to process'

    os.chdir('..')
