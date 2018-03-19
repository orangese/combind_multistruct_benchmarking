import os
import sys

from parse_chembl import load_chembl_proc, load_drugs

from schrodinger.structure import StructureReader, StructureWriter

def all_mcss(num_chembl=0):

    drugs = load_drugs()

    chembl_info = {'{}_lig'.format(l):lig for l, lig in load_chembl_proc().items()}
    all_ligs = sorted([l.split('.')[0] for l in os.listdir('ligands/unique')])
    pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL'] + [l for l in all_ligs if l.split('_')[0] in drugs]
    chembl_ligs = [l for l in all_ligs if l in chembl_info and l not in drugs and chembl_info[l].valid_stereo]
    chembl_ligs.sort(key=lambda x: chembl_info[x].ki)

    #print pdb_ligs
    #print chembl_ligs[:num_chembl]
    #return

    os.system('mkdir -p ligands/mcss')

    os.system('rm -f ligands/mcss/*_in.mae')

    all_pairs = []

    for i1 in range(min(50,len(pdb_ligs))):
        for i2 in range(i1+1,min(50,len(pdb_ligs))):
            l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
            if not os.path.exists('ligands/mcss/{}-{}.mae'.format(l1, l2)):
                all_pairs.append((l1, l2))

    for l1 in pdb_ligs:
        for l2 in chembl_ligs[:num_chembl]:
            if not os.path.exists('ligands/mcss/{}-{}.mae'.format(l1, l2)):
                all_pairs.append((l1, l2))

    for i1 in range(min(num_chembl,len(chembl_ligs))):
        for i2 in range(i1+1,min(num_chembl,len(chembl_ligs))):
            l1, l2 = chembl_ligs[i1], chembl_ligs[i2]
            if not os.path.exists('ligands/mcss/{}-{}.mae'.format(l1, l2)):
                all_pairs.append((l1, l2))

    if len(all_pairs) == 0: return

    if len(all_pairs) > 0:#5000:
        print len(all_pairs), 'mcss pairs left'

    with open('ligands/mcss/mcss_in.sh', 'w') as f:
        f.write('#!/bin/bash\nmodule load schrodinger\n')

        for (l1, l2) in all_pairs[:5000]:
            name = '{}-{}'.format(l1, l2) 
            st1 = StructureReader('ligands/unique/{}.mae'.format(l1)).next()
            st2 = StructureReader('ligands/unique/{}.mae'.format(l2)).next()
            st1._setTitle(l1)
            st2._setTitle(l2)

            st_wr = StructureWriter('ligands/mcss/{}_in.mae'.format(name))
            st_wr.append(st1)
            st_wr.append(st2)
            st_wr.close()
                    
            f.write('$SCHRODINGER/utilities/canvasMCS -JOB {} -WAIT -imae {}_in.mae -omae {}.mae -atomtype 13 -nobreakring -exclusive\n'.format(name, name, name))

    os.chdir('ligands/mcss')
    os.system('sbatch -p rondror -n 10 -t 1:00:00 -o mcss.out mcss_in.sh')
    os.chdir('../..')

