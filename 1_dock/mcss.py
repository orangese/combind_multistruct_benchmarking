import os
import sys
import itertools

from parse_chembl import load_chembl_proc, load_drugs

from schrodinger.structure import StructureReader, StructureWriter

atom_typing = 7
out_dir = 'ligands/mcss/mcss{}'.format(atom_typing)

group_size=180
def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

def all_mcss(num_chembl=0):

    drugs = load_drugs()

    chembl_info = {'{}_lig'.format(l):lig for l, lig in load_chembl_proc().items()}
    all_ligs = sorted([l.split('.')[0] for l in os.listdir('ligands/unique')])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL'] + [l for l in all_ligs if l.split('_')[0] in drugs]
    chembl_ligs = [l for l in all_ligs if l in chembl_info and l not in drugs]# and chembl_info[l].valid_stereo]
    chembl_ligs.sort(key=lambda x: chembl_info[x].ki)

    #print pdb_ligs
    #print chembl_ligs[:num_chembl]
    #return

    os.system('mkdir -p {}'.format(out_dir))

    #os.system('rm -f {}/*_in.mae'.format(out_dir))
    for f in os.listdir(out_dir):
        if f.split('_')[-1] == 'in.mae':
            os.system('rm -f {}/{}'.format(out_dir, f))

    all_pairs = []

    for i1 in range(min(50,len(pdb_ligs))):
        for i2 in range(i1+1,min(50,len(pdb_ligs))):
            l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
            if not os.path.exists('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                all_pairs.append((l1, l2))

    for l1 in pdb_ligs:
        for l2 in chembl_ligs[:num_chembl]:
            if not os.path.exists('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                all_pairs.append((l1, l2))

    for i1 in range(min(num_chembl,len(chembl_ligs))):
        for i2 in range(i1+1,min(num_chembl,len(chembl_ligs))):
            l1, l2 = chembl_ligs[i1], chembl_ligs[i2]
            if not os.path.exists('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                all_pairs.append((l1, l2))

    if len(all_pairs) == 0: return

    if len(all_pairs) > 0:#5000:
        print len(all_pairs), 'mcss pairs left'

    os.chdir(out_dir)
    os.system('rm -f *.sh')
    for i, pairs in enumerate(grouper(group_size, all_pairs)):
        with open('{}mcss.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2) 
                st1 = StructureReader('../../unique/{}.mae'.format(l1)).next()
                st2 = StructureReader('../../unique/{}.mae'.format(l2)).next()
                st1._setTitle(l1)
                st2._setTitle(l2)

                st_wr = StructureWriter('{}_in.mae'.format(name))
                st_wr.append(st1)
                st_wr.append(st2)
                st_wr.close()
                    
                f.write('$SCHRODINGER/utilities/canvasMCS -JOB {} -WAIT -imae {}_in.mae -omae {}.mae -atomtype {} -nobreakaring -exclusive\n'.format(name, name, name, atom_typing))
            f.write('wait\n')
        os.system('sbatch -p owners --tasks=6 --cpus-per-task=1 -t 1:30:00 -o mcss.out {}mcss.sh'.format(i))
    os.chdir('../../..')

