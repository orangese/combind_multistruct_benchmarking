import os
import sys
import itertools

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty

atom_typing = 7
out_dir = 'ligands/mcss/mcss{}'.format(atom_typing)

command = '$SCHRODINGER/utilities/canvasMCS -JOB {} -WAIT -imae {}_in.mae -omae {}.mae -atomtype {} -nobreakaring -exclusive -stop 10'
command_csv = '$SCHRODINGER/utilities/canvasMCS -JOB {} -WAIT -imae {}_in.mae -ocsv {}.csv -atomtype {} -nobreakaring -exclusive -stop 10'

l_dir = 'ligands/prepared_ligands'
l_path = '{}/{}/{}.mae'

queue = 'rondror'

group_size=15
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

def init_mcss(lm, chembl=None, max_num=20):
    os.system('mkdir -p {}'.format(out_dir))

    all_pairs = set([])
    if chembl is not None:
        for f, f_data in chembl.items():
            for q,c in f_data.items():
                for i in range(len(c)):
                    for j in range(i+1,len(c)):
                        all_pairs.add((c[i],c[j]))
    else:
        pdb_ligs = lm.pdb[:max_num]
        for i1 in range(len(pdb_ligs)):
            for i2 in range(i1+1,len(pdb_ligs)):
                l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
                all_pairs.add((l1, l2))

        for l1 in pdb_ligs:
            for l2 in lm.chembl():
                #txt_out = '{}/{}-{}.txt'.format(out_dir, l1, l2)
                #if os.path.exists(txt_out) and os.stat(txt_out).st_size == 0:
                #    os.system('rm {}'.format(txt_out))
                #if not os.path.exists(txt_out):
                all_pairs.add((l1, l2))

    mcss_compute(all_pairs)

def mcss_compute(all_pairs):
    no_mcss = []
    no_size_file = []

    for l1,l2 in all_pairs:

        csv_out = '{}/{}-{}.csv'.format(out_dir, l1, l2)
        mae_out = '{}/{}-{}.mae'.format(out_dir, l1, l2)
        txt_out = '{}/{}-{}.txt'.format(out_dir, l1, l2)

        if os.path.exists(txt_out):

            if os.stat(txt_out).st_size == 0:
                print l1,l2, 'empty file found'
                os.system('rm {}'.format(txt_out))

            elif os.path.exists(mae_out) or os.path.exists(csv_out):
                # delete the structure file once the text file is written
                os.system('rm -f {} {}'.format(mae_out, csv_out))
                os.system('rm -f {}/{}-{}_in.mae'.format(out_dir, l1, l2))

        elif not os.path.exists(mae_out) and not os.path.exists(csv_out):
            no_mcss.append((l1,l2))
        elif not os.path.exists('{}/{}-{}.txt'.format(out_dir, l1, l2)):
            no_size_file.append((l1,l2))

    if len(no_mcss) > 0:
        print len(no_mcss), 'mcss pairs left'
        get_mcss(no_mcss)
    if len(no_size_file) > 0:
        print len(no_size_file), 'size pairs left'
        get_size_file(no_size_file)

def get_mcss(mcss_pairs):
    for i, pairs in enumerate(grouper(group_size, mcss_pairs)):
        with open('{}/cmcss{}.sh'.format(out_dir,i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2) 
                st1 = StructureReader(l_path.format(l_dir,l1,l1)).next()
                st2 = StructureReader(l_path.format(l_dir,l2,l2)).next()
                st1.title = l1
                st2.title = l2

                st_wr = StructureWriter('{}/{}_in.mae'.format(out_dir,name))
                st_wr.append(st1)
                st_wr.append(st2)
                st_wr.close()
                    
                f.write(command_csv.format(name, name, name, atom_typing)+'\n')
            f.write('wait\n')
        os.chdir(out_dir)
        os.system('sbatch -p {} --tasks=2 --cpus-per-task=1 --ntasks-per-socket=2 -t 1:00:00 -o mcss.out cmcss{}.sh'.format(queue,i))
        os.chdir('../../..')

def get_size_file(size_pairs):
    for l1,l2 in size_pairs:
        csv_out = '{}/{}-{}.csv'.format(out_dir, l1, l2)
        mae_out = '{}/{}-{}.mae'.format(out_dir, l1, l2)

        try:
            with open('{}/{}-{}.txt'.format(out_dir, l1, l2),'w') as f:
                if os.path.exists(mae_out):
                    for st in StructureReader('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                        prop = _StructureProperty(st)
                        smarts = prop['s_canvas_MCS_SMARTS']
                        msize = prop['i_canvas_MCS_Atom_Count'] 
                        size = len([a for a in st.atom if a.element != 'H'])
                    f.write('{},{},{},{}\n'.format(st.title, size, msize, smarts))
                elif os.path.exists(csv_out):
                    for i,st in enumerate(StructureReader('{}/{}-{}_in.mae'.format(out_dir, l1, l2))):
                        size = len([a for a in st.atom if a.element != 'H'])
                        with open(csv_out) as csv_f:
                            for j,line in enumerate(csv_f):
                                if i+1 == j:
                                    line = line.strip().split(',')
                                    lig = line[1]
                                    msize = int(line[5])
                                    smarts = line[-1]
                                    assert lig == st.title
                                    f.write('{},{},{},{}\n'.format(st.title, size, msize, smarts))
                                    #print '{},{},{},{}\n'.format(st.title, size, msize, smarts)

        except Exception as e:
            print e
            print 'size file error', l1, l2
            #os.system('rm {}/{}-{}.mae'.format(out_dir, l1, l2))
            os.system('rm {}/{}-{}.txt'.format(out_dir, l1, l2))
        #break


