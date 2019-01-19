import os
import sys

def pick_helpers(lm, maxnum=21):

    parent = 'chembl/helpers'
    os.system('mkdir -p {}'.format(parent))

    stereo_filter = lambda c,ci: ci[c].valid_stereo

    ki_sort = lambda c: lm.chembl_info[c].ki

    # filters
    all_options = [
        'best_affinity.txt',
        'best_mcss.txt',
    ]

    num_chembl = 30
    for f in all_options:
        fpath = '{}/{}'.format(parent, f)
        if not os.path.exists(fpath):
            print('picking chembl ligands', f)
            chembl_ligs = lm.chembl()
            with open(fpath,'w') as fi:
                for q in lm.get_xdocked_ligands(maxnum):
                    # sort and remove duplicates
                    print(q)
                    if 'mcss' in f:
                        lm.mcss.load_mcss()
                        sorted_helpers = sorted(chembl_ligs, key=ki_sort)
                        sorted_helpers = lm.mcss.sort_by_mcss(q, sorted_helpers)
                    elif 'affinity' in f:
                        sorted_helpers = sorted(chembl_ligs, key=ki_sort)
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
                q, chembl = line.strip().split(':')
                helpers[fname][q] = chembl.split(',')
    return helpers
