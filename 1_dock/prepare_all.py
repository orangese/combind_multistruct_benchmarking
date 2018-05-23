import os
import sys

from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files
from init_mcss import init_mcss

from grids import make_grids
from dock import dock
from fp import fp
from mcss import mcss

from chembl_sort import get_ligands, proc_ligands
from pick_helpers import pick_helpers, load_helpers

os.chdir('../../data')

datasets = sys.argv[1:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir('.')) if d[0] != '.' and d[-3:] != 'old']

#grids = {'D2R':'6CM4','AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
#         'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23','TRPV1':'3J5Q','SIGMA1':'5HK1'}

grids = {
    'D2R':'6CM4',
    'AR':'2PNU',
    'B1AR':'2VT4',
    'TRPV1':'3J5Q',
    'SIGMA1':'5HK1',
    '5HT2B':'4IB4',
    'DTRANSP':'4M48', 
    'M3':'4U15'
}

for i, d in enumerate(datasets):

    print i, d

    os.chdir(d)

    pdb_st = None
    if os.path.exists('structures/downloads'):    
        pdb_st = sort_downloads()
     
    # 1. prepare proteins    
    process_structs()
    align_structs()
    sort_files()
    make_grids()
    dock(os.listdir('docking/grids'))
     
    # 2. prepare ligands
    get_ligands()    
    proc_ligands()
    init_mcss()

    if False: # force redo of chembl info (do this if new chembl ligands have been added)
        os.system('rm -f chembl/helpers/*')
        os.system('rm -f chembl/duplicates.txt')
        os.system('rm -f chembl/molw.txt')

    # 3. decide what ligands to use
    pick_helpers()
    h = load_helpers()
    init_mcss(h)

    # 4. dock/fp/mcss those ligands
    if len(os.listdir('docking/grids')) == 0: 
        os.chdir('..')
        continue

    to_prep = [sorted(os.listdir('docking/grids'))[0]]
    if d in grids and grids[d] not in to_prep:
        to_prep += grids[d]

    dock(to_prep, h)
    fp(to_prep)
    mcss(to_prep, h)    


    os.chdir('..')











