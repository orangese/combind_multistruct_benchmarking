import os
from multiprocessing import Pool

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def processHelper(struct):
    options = "-WAIT -fillsidechains -f 3 -fix -samplewater -delwater_hbond_cutoff 2 -keepfarwat -captermini"  
    command = "{}/utilities/prepwizard {} -j temp-{} ../stripped/{}.mae {}.mae".format(SCHRODINGER, options, struct, struct, struct)
    os.system(command)

def process():   
    os.system('mkdir -p processed')
    os.chdir('processed')
    
    failed = []
    for struct in [f.split('.')[0] for f in os.listdir('../stripped')]:
        for i in range(5):
            if '{}.mae'.format(struct) in os.listdir('.'): 
                print '{} succeeded!'.format(struct)
                break 
            print 'attempting to process {}. attempt {}'.format(struct, i+1)
            processHelper(struct)
        else:
            failed.append(struct)

    print '{} failed files'.format(len(failed))
    print failed

    os.system('rm temp* *missing* *.out')
    os.chdir('..')
