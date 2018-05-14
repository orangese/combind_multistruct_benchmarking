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
    if not os.path.exists('ligands/unique'): return
    from schrodinger.structure import StructureReader
    with open('chembl/molw.txt', 'w') as fname:
        for f in os.listdir('ligands/unique'):
            st = StructureReader('ligands/unique/{}'.format(f)).next()
            mw = st.total_weight
            fname.write('{},{}\n'.format(f.split('.')[0],mw))

def write_duplicates():
    if not os.path.exists('ligands/unique'): return
    duplicates = []
    from schrodinger.structure import StructureReader
    all_ligs = {f.split('.')[0] : StructureReader('ligands/unique/'+f).next() for f in os.listdir('ligands/unique')}
    for l1 in [l for l in all_ligs if l[:6] != 'CHEMBL']:
        for l2 in [l for l in all_ligs if l[:6] == 'CHEMBL']:
            if all_ligs[l1].isEquivalent(all_ligs[l2]):
                duplicates.append((l1, l2))
    with open('chembl/duplicates.txt','w') as f:
        for l1,l2 in duplicates:
            f.write('{},{}\n'.format(l1,l2))

def read_duplicates(dir_path=None):
    duplicates = []
    fpath = 'chembl/duplicates.txt'
    if dir_path is not None:
        fpath = '{}/{}'.format(dir_path, fpath)
    if not os.path.exists(fpath) and os.path.exists('chembl'):
        print 'writing duplicates'
        write_duplicates()
    if not os.path.exists(fpath): return {}
    with open(fpath) as f:
        for line in f:
            l1, l2 = line.strip().split(',')
            duplicates.append((l1,l2))
    return duplicates

def read_mcss(dir_path=None):
    mcss = {}

    ligpath = 'ligands/unique'
    mcsspath = 'ligands/mcss/mcss7'
    if dir_path is not None:
        ligpath = '{}/{}'.format(dir_path, ligpath)
        mcsspath = '{}/{}'.format(dir_path, mcsspath)

    if not os.path.exists(ligpath) or not os.path.exists(mcsspath): 
        print os.listdir('.')
        return {} 

    all_ligs = sorted([l.split('.')[0] for l in os.listdir(ligpath)])
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







