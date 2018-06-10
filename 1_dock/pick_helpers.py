import os
import sys

from chembl_props import read_mcss

def pick_helpers(lm, maxnum=20):

    parent = 'chembl/helpers'
    os.system('mkdir -p {}'.format(parent))

    stereo_filter = lambda c,ci: ci[c].valid_stereo

    mcss_sort = lambda c: -mcss[q].get(c,(0,0,0))[2]
    ki_sort = lambda c: lm.chembl_info[c].ki

    # sorting function plus (optional) additional filters
    all_options = {
        'best_affinity.txt': [ki_sort],
        'best_mcss.txt': [mcss_sort],
        'best_affinity_plus_stereo.txt': [ki_sort, stereo_filter],
        'best_mcss_plus_stereo.txt': [mcss_sort, stereo_filter]
    }

    queries = lm.docked(lm.pdb)[:maxnum]
    mcss = { q : read_mcss(containing_lig=q) for q in queries }

    num_chembl = 30    
    for f, conditions in all_options.items():
        fpath = '{}/{}'.format(parent, f)
        if not os.path.exists(fpath):
            print 'picking chembl ligands', f
            # apply filters
            chembl_ligs = lm.chembl(conditions[1:])
            with open(fpath,'w') as fi:
                for q in queries:
                    # sort and remove duplicates
                    unique = lm.unique([q]+sorted(chembl_ligs,key=conditions[0]))
                    fi.write('{}:{}\n'.format(q, ','.join(unique[1:num_chembl+1])))

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



