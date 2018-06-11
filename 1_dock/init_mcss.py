import os
import sys
from grouper import grouper

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty

command = '$SCHRODINGER/utilities/canvasMCS -imae {}_in.mae -omae {}.mae -atomtype 7 -exclusive -stop 10'
command_csv = '$SCHRODINGER/utilities/canvasMCS -imae {}_in.mae -ocsv {}.csv -atomtype 7 -exclusive -stop 10'

l_path = 'ligands/prepared_ligands/{}/{}_neutral.mae'

mae_in  = 'ligands/mcss/{}/{}-{}_in.mae'
csv_out = 'ligands/mcss/{}/{}-{}.csv'
mae_out = 'ligands/mcss/{}/{}-{}.mae'
txt_out = 'ligands/mcss/{}/{}-{}.txt'

queue = 'rondror'
group_size=10

def init_mcss(lm, chembl=None, max_num=20):
    os.system('mkdir -p ligands/mcss/{}'.format(lm.sp['mcss']))

    all_pairs = set([])
    if chembl is not None:
        for f, f_data in chembl.items():
            for q,c in f_data.items():
                for i in range(len(c)):
                    for j in range(i+1,len(c)):
                        all_pairs.add((c[i],c[j]))
    else:
        pdb_ligs = [l for l in lm.pdb if os.path.exists(l_path.format(l,l))][:max_num]
        for i1 in range(len(pdb_ligs)):
            for i2 in range(i1+1,len(pdb_ligs)):
                l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
                all_pairs.add((l1, l2))

        for l1 in pdb_ligs:
            for l2 in lm.chembl():
                all_pairs.add((l1, l2))

    mcss_compute(all_pairs, lm.sp['mcss'])

def mcss_compute(all_pairs, mdir):

    no_mcss = []
    no_size_file = []

    for l1,l2 in all_pairs:

        if os.path.exists(txt_out.format(mdir,l1,l2)):
            if os.path.exists(mae_in.format(mdir,l1,l2)):
                print mae_in.format(mdir,l1,l2), 'error, retry'
                #os.system('rm -f {}/{}-{}*'.format(out_dir, l1, l2))
        else:
            if not os.path.exists(mae_out.format(mdir,l1,l2)) and not os.path.exists(csv_out.format(mdir,l1,l2)):
                no_mcss.append((l1,l2))
            else:
                no_size_file.append((l1,l2))

    if len(no_mcss) > 0:
        print len(no_mcss), 'init_mcss pairs left, step1'
        get_mcss(no_mcss, mdir)
    if len(no_size_file) > 0:
        print len(no_size_file), 'init_mcss pairs left, step2'
        get_size_file(no_size_file, mdir)

def get_mcss(mcss_pairs,mdir):
    for i, pairs in enumerate(grouper(group_size, mcss_pairs)):
        with open('ligands/mcss/{}/cmcss{}.sh'.format(mdir,i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2) 
                st1 = StructureReader(l_path.format(l1,l1)).next()
                st2 = StructureReader(l_path.format(l2,l2)).next()
                st1.title = l1
                st2.title = l2

                st_wr = StructureWriter(mae_in.format(mdir,l1,l2))
                st_wr.append(st1)
                st_wr.append(st2)
                st_wr.close()
                    
                f.write(command_csv.format(name, name, name)+'\n')
            f.write('wait\n')
        os.chdir('ligands/mcss/{}'.format(mdir))
        #os.system('sbatch -p {} --tasks=1 --cpus-per-task=1 -t 1:00:00 -o mcss.out cmcss{}.sh'.format(queue,i))
        os.chdir('../../..')

def get_size_file(size_pairs,mdir):
    for l1,l2 in size_pairs:
        try:
            line_count = 0
            with open(txt_out.format(mdir,l1,l2),'w') as f:
                if os.path.exists(mae_out.format(mdir,l1,l2)):
                    for st in StructureReader(mae_out.format(mdir,l1,l2)):
                        prop = _StructureProperty(st)
                        smarts = prop['s_canvas_MCS_SMARTS']
                        msize = prop['i_canvas_MCS_Atom_Count'] 
                        size = len([a for a in st.atom if a.element != 'H'])
                    f.write('{},{},{},{}\n'.format(st.title, size, msize, smarts))
                    line_count += 1
                elif os.path.exists(csv_out.format(mdir,l1,l2)):
                    for i,st in enumerate(StructureReader(mae_in.format(mdir,l1,l2))):
                        size = len([a for a in st.atom if a.element != 'H'])
                        with open(csv_out.format(mdir,l1,l2)) as csv_f:
                            for j,line in enumerate(csv_f):
                                if i+1 == j:
                                    line = line.strip().split(',')
                                    lig = line[1]
                                    msize = int(line[5])
                                    smarts = line[-1]
                                    assert lig == st.title
                                    f.write('{},{},{},{}\n'.format(st.title, size, msize, smarts))
                                    line_count += 1
            assert line_count == 2, mae_in.format(mdir,l1,l2)
            os.system('rm -f {} {} {}'.format(csv_out.format(mdir,l1,l2), mae_out.format(mdir,l1,l2), mae_in.format(mdir,l1,l2)))
        except Exception as e:
            print e
            print 'size file error', l1, l2, line_count
            os.system('rm -f ligands/mcss/{}/{}-{}.txt'.format(mdir, l1, l2))


