import os
import sys

def pick_helpers(lm, maxnum=20):

    parent = 'chembl/helpers'
    os.system('mkdir -p {}'.format(parent))

    stereo_filter = lambda c,ci: ci[c].valid_stereo

    ki_sort = lambda c: lm.chembl_info[c].ki

    # filters
    all_options = {
        'best_affinity.txt': [],#ki_sort],
        'best_mcss.txt': []#mcss_sort]
        #'best_affinity_plus_stereo.txt': [ki_sort, stereo_filter],
        #'best_mcss_plus_stereo.txt': [mcss_sort, stereo_filter]
    }

    num_chembl = 30    
    for f, filters in all_options.items():
        fpath = '{}/{}'.format(parent, f)
        if not os.path.exists(fpath):
            print 'picking chembl ligands', f
            # apply filters
            chembl_ligs = lm.chembl(filters)
            with open(fpath,'w') as fi:
                for q in lm.docked(lm.pdb)[:maxnum]:
                    # sort and remove duplicates
                    if 'mcss' in f:
                        lm.mcss.load_mcss(set([q]), set(chembl_ligs))
                        sorted_helpers = lm.mcss.sort_by_mcss(q, set(chembl_ligs),lambda x: x.m_sz)
                    elif 'affinity' in f:
                        sorted_helpers = sorted(chembl_ligs,key=ki_sort)
                    unique = lm.unique([q]+sorted_helpers)
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



