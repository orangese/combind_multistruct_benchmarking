import os
import sys

from parse_chembl import load_chembl_proc


def pick_helpers():
    chembl_info = load_chembl_proc()
    all_ligs = sorted([l.split('.')[0] for l in os.listdir('ligands/unique')])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL']
    # most basic filtering: ki < 1 uM and molecular weight < 1000 Da
    chembl_ligs = [l for l in all_ligs if l in chembl_info and chembl_info[l].ki <= 1000 and chembl_info[l].mw <= 1000]

    num_chembl = 30

    parent = 'chembl/helpers'
    os.system('mkdir -p {}'.format(parent))

    all_options = {
        'best_affinity.txt': lambda c: chembl_info[c].ki,
        'best_mcss.txt': lambda c: -chembl_info[c].mcss.get(q,(0,0,0))[2]
    }
    
    for f, sort_f in all_options.items():
        fpath = '{}/{}'.format(parent, f)
        if not os.path.exists(fpath):
            with open(fpath,'w') as fi:
                for q in pdb_ligs:
                    chembl_ligs.sort(key=sort_f)
                    fi.write('{}:{}\n'.format(q, ','.join(chembl_ligs[:num_chembl])))

def load_helpers(fname, dirpath=None):
    fpath = 'chembl/helpers/{}'.format(fname)
    if dirpath is not None:
        fpath = '{}/{}'.format(dirpath, fpath)

    helpers = {}
    with open(fpath) as f:
        for line in f:
            q,chembl = line.strip().split(':')
            helpers[q] = chembl.split(',')

    return helpers



