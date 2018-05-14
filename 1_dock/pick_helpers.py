import os
import sys

from parse_chembl import load_chembl_proc


def pick_helpers():

    parent = 'chembl/helpers'
    os.system('mkdir -p {}'.format(parent))

    stereo_filter = lambda c: chembl_info[c].valid_stereo
    mw_filter = lambda c: chembl_info[c].mw <= 1000
    ki_filter = lambda c: chembl_info[c].ki <= 1000 

    mcss_sort = lambda c: -chembl_info[c].mcss.get(q,(0,0,0))[2]
    ki_sort = lambda c: chembl_info[c].ki

    # sorting function plus (optional) additional filters
    all_options = {
        'best_affinity.txt': [ki_sort],
        'best_mcss.txt': [mcss_sort],
        'best_affinity_plus_stereo.txt': [ki_sort, stereo_filter],
        'best_mcss_plus_stereo.txt': [mcss_sort, stereo_filter]
    }

    for f in all_options:
        if not os.path.exists('{}/{}'.format(parent, f)):
            os.system('rm -f chembl/helpers/*')
            break
    else: return

    chembl_info = load_chembl_proc(load_mcss=True)
    all_ligs = sorted([l.split('.')[0] for l in os.listdir('ligands/unique')])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL']

    num_chembl = 30
    
    for f, conditions in all_options.items():
        fpath = '{}/{}'.format(parent, f)
        if not os.path.exists(fpath):
            filters = [mw_filter, ki_filter] # apply to all
            if len(conditions) > 1: filters.extend(conditions[1:])
            chembl_ligs = [l for l in all_ligs if l in chembl_info and False not in [f(l) for f in filters]]
            with open(fpath,'w') as fi:
                for q in pdb_ligs:
                    chembl_ligs.sort(key=conditions[0])
                    fi.write('{}:{}\n'.format(q, ','.join(chembl_ligs[:num_chembl])))

def load_helpers(dirpath=None):
    fpath = 'chembl/helpers'
    if dirpath is not None:
        fpath = '{}/{}'.format(dirpath, fpath)

    helpers = {}
    for fname in os.listdir(fpath):
        if fname[0] == '.' or fname.split('.')[-1] != 'txt': continue
        helpers[fname] = {}
        with open('{}/{}'.format(fpath, fname)) as f:
            for line in f:
                q,chembl = line.strip().split(':')
                helpers[fname][q] = chembl.split(',')

    return helpers



