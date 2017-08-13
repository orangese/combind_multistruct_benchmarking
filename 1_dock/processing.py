import os
from multiprocessing import Pool
import slurm

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(struct):
    if '{}.mae'.format(struct) in os.listdir('.'):
        return struct, True
    options = "-WAIT -fillsidechains -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini"  
    command = "{}/utilities/prepwizard {} -j temp-{} ../stripped/{}.mae {}.mae".format(SCHRODINGER, options, struct, struct, struct)
    slurm.salloc(command, '1','1:00:00')
    return struct, '{}.mae'.format(struct) in os.listdir('.')

def process():   
    os.system('mkdir -p processed')
    os.chdir('processed')
    
    unfinished_structs = [f.split('.')[0] for f in os.listdir('../stripped')]
    pool = Pool(len(unfinished_structs))

    for i in range(5):
        processing_structs = unfinished_structs
        unfinished_structs = []
        print 'iteration {}, processing {} structs'.format(i+1, len(processing_structs))
        print processing_structs
        for s, done in pool.imap_unordered(processHelper, processing_structs):
            if done: print s, ' succeeded!'
            else: unfinished_structs.append(s)
        if len(unfinished_structs) == 0: 
            print 'all done!'
            break
    else:
        print unfinished_structs, 'failed to process'

    os.system('rm temp* *missing* *.out')
    os.chdir('..')
