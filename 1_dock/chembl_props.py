import os
import sys

def read_molw(dir_path=None):
    molw = {}
    fpath = 'chembl/molw.txt'
    if dir_path is not None:
        fpath = '{}/{}'.format(dir_path, fpath)
    if not os.path.exists(fpath) and os.path.exists('chembl'):
        print 'writing molw'
        write_molw()
    if not os.path.exists(fpath): return {}
    with open(fpath) as fname:
        for line in fname:
            chembl, mw = line.strip().split(',')
            molw[chembl] = float(mw)
    return molw

def write_molw():
    l_dir = 'ligands/prepared_ligands'
    l_path = '{}/{}/{}.mae'

    if not os.path.exists(l_dir): return
    from schrodinger.structure import StructureReader
    with open('chembl/molw.txt', 'w') as fname:
        for f in os.listdir(l_dir):
            if not os.path.exists(l_path.format(l_dir, f, f)): continue
            st = StructureReader(l_path.format(l_dir, f,f)).next()
            mw = st.total_weight
            fname.write('{},{}\n'.format(f,mw))

def write_duplicates():
    l_dir = 'ligands/prepared_ligands'
    l_path = '{}/{}/{}.mae'

    if not os.path.exists(l_dir): return
    duplicates = {}
    from schrodinger.structure import StructureReader
    all_ligs = [l for l in os.listdir(l_dir) if os.path.exists(l_path.format(l_dir, l, l))]
    all_ligs = {l : StructureReader(l_path.format(l_dir,l,l)).next() for l in all_ligs}
    for l1 in all_ligs:#[l for l in all_ligs if l[:6] != 'CHEMBL']:
        for l2 in duplicates:
            if all_ligs[l1].isEquivalent(all_ligs[l2], True):
                duplicates[l2].append(l1)
                break
        else:
            duplicates[l1] = [l1]
    with open('chembl/duplicates.txt','w') as f:
        for l_list in duplicates:
            f.write('{}\n'.format(','.join(l_list)))

def read_duplicates(dir_path=None):
    duplicates = {}
    fpath = 'chembl/duplicates.txt'
    if dir_path is not None:
        fpath = '{}/{}'.format(dir_path, fpath)
    if not os.path.exists(fpath) and os.path.exists('chembl'):
        print 'writing duplicates'
        write_duplicates()
    if not os.path.exists(fpath): return {}
    with open(fpath) as f:
        for line in f:
            l_list = line.strip().split(',')
            duplicates[l_list[0]] = l_list#.append((l1,l2))
    return duplicates

def read_mcss(dir_path=None):
    mcss = {}

    l_dir = 'ligands/prepared_ligands'
    l_path = '{}/{}/{}.mae'
    mcsspath = 'ligands/mcss/mcss7'
    if dir_path is not None:
        l_dir = '{}/{}'.format(dir_path, l_dir)
        mcsspath = '{}/{}'.format(dir_path, mcsspath)

    if not os.path.exists(l_dir) or not os.path.exists(mcsspath): 
        print os.listdir('.')
        return {} 

    all_ligs = sorted([l for l in os.listdir(l_dir) if os.path.exists(l_path.format(l_dir, l, l))])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL']
    for q in pdb_ligs:
        qpairs = [f for f in os.listdir(mcsspath) if f.split('-')[0] == q and f[-3:] == 'txt']
        for pfile in qpairs:
            try:
                with open('{}/{}'.format(mcsspath,pfile)) as f:
                    l1,s1,m1,smarts1 = f.readline().strip().split(',')
                    l2,s2,m2,smarts2 = f.readline().strip().split(',')
                    assert l1 == q, '{} {} {} error'.format(l1,q,pfile)
                    assert m1 == m2, 'msize error {}'.format(pfile)
                    if l1 not in mcss: mcss[l1] = {}
                    mcss[l1][l2] = (int(s1),int(s2),int(m1))
            except Exception as e:
                print e
                print pfile, l1,l2
                os.system('rm -f {}/{}'.format(mcsspath,pfile))

    return mcss







