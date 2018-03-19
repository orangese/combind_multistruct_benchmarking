import os
import sys

from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files

from unique_ligands import filter_duplicates
from ligand_features import output_features
from mcss import all_mcss

from core_proc import core_proc
from grids import make_grids
from dock import dock_dataset, core_dock

from chembl_sort import get_ligands, get_drugs, proc_ligands

os.chdir('../../data')

datasets = sys.argv[1:]
if len(datasets) == 0: datasets = sorted(os.listdir('.'))

grids = {'D2R':'6CM4','AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
         'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23'}

for i, d in enumerate(datasets):

    #if d[0] == '.': continue

    if d not in grids: continue

    # INPUT:
    # <protein>/structures/raw_files
    # should contain <lig>_lig.mae and <prot>_prot.mae
    # <protein>/chembl/<protein>_<organism>.xls

    print i, d

    os.chdir(d)

    pdb_st = None
    if os.path.exists('structures/downloads'):    
        pdb_st = sort_downloads()
    
    process_structs()
    align_structs()#True)
    sort_files()
    make_grids()
    #continue
    # process ligands
    get_ligands()
    if os.path.exists('chembl/drugs.txt'):
        get_drugs()

    proc_ligands()
    filter_duplicates()
    output_features()
    all_mcss(num_chembl=25)

    #continue

    #grids = sorted([st for st in os.listdir('structures/aligned_files')
    #                if st[0] != '.' and st+'_lig.mae' not in os.listdir('ligands/duplicate')])

    #if os.path.exists('docking/core'):
    #    core_proc()
    #    core_dock(grids)

    dock_dataset([grids[d]], num_chembl=25)
    
    os.chdir('..')


