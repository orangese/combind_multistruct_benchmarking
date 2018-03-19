import os
import sys

import itertools

data = '/scratch/PI/rondror/jbelk/method/data'

os.chdir(data)

glide_dir = sys.argv[1]
out_dir = sys.argv[2]

d_list = sys.argv[3:]

if len(d_list) == 0:
    d_list = sorted(os.listdir('.'))

datasets = {'D2R':'6CM4','AR':'2PNU','A2AR':'2YDO','B1AR':'2VT4','B2AR':'2RH1','CHK1':'2BRN', 'PLK1':'2OWB',
            'VITD':'2HB7','BRAF':'3IDP','JAK2':'3KRR','CDK2':'1H1S','ERA':'1A52','GCR':'3K23'}

script = '/scratch/PI/rondror/jbelk/method/code/2_ifp/proc_mcss.py'

def grouper(n, iterable, fillvalue=None):
    #"grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue) #NOTE: izip_longest is zip_longest in python2

for i, d in enumerate(d_list):
    if d not in datasets: continue

    print i, d

    st = datasets[d]

    os.chdir(d)

    os.system('mkdir -p mcss/{}/{}'.format(out_dir, st))
    os.system('rm -f mcss/{}/{}/*.sh'.format(out_dir, st))
    for fname in os.listdir('mcss/{}/{}'.format(out_dir, st)):
        if os.path.getsize('mcss/{}/{}/{}'.format(out_dir, st, fname)) == 0:
            print 'found empty file', fname
            os.system('rm -f mcss/{}/{}/{}'.format(out_dir, st, fname))

    for fname in os.listdir('mcss/{}/{}'.format(out_dir, st)):
        if fname[-3:] == 'out':
            with open('mcss/{}/{}/{}'.format(out_dir, st, fname)) as f:
                for line in f:
                    line = line.strip().split(' ')
                    if line[0] == 'error?':
                        print line[1]
                        #print f.read()
                        #os.system('rm mcss/{}/{}/{}*'.format(out_dir, st, line[1]))

    #os.chdir('..')
    #continue    

    all_mcss = sorted([m.split('.')[0] for m in os.listdir('ligands/mcss') if m.split('.')[-1] == 'mae'])
    all_mcss = [tuple(m.split('-')) for m in all_mcss if m[-3:] != '_in']

    pdb_ids = [l.split('.')[0] for l in os.listdir('ligands/unique')]
    pdb_ids = sorted([l for l in pdb_ids if l+'.mae' in os.listdir('structures/ligands')])
    pdb_ids = set(pdb_ids)

    unfinished_pairs = []
    for l1, l2 in all_mcss:
        #print l1,l2
        pv1, pv2 = '{}-to-{}'.format(l1, st), '{}-to-{}'.format(l2, st)
        if not os.path.exists('docking/{}/{}/{}_pv.maegz'.format(glide_dir, pv1, pv1)): continue
        if not os.path.exists('docking/{}/{}/{}_pv.maegz'.format(glide_dir, pv2, pv2)): continue
        if os.path.exists('mcss/{}/{}/{}-{}.csv'.format(out_dir, st, l1, l2)): continue
        #if l1 not in pdb_ids or l2 not in pdb_ids: continue
        #if l1[:6] == 'CHEMBL' or l2[:6] == 'CHEMBL': continue
        #print l1,l2
        unfinished_pairs.append((l1, l2))

    if len(unfinished_pairs) > 0:
        print len(unfinished_pairs)#, unfinished_pairs
    #break

    os.chdir('mcss/{}/{}'.format(out_dir, st))
    os.system('rm -f *.in')
    for i, pairs in enumerate(grouper(4, unfinished_pairs)):
    #for l1, l2 in unfinished_pairs:
        #print l1, l2
        with open('{}_in.sh'.format(i), 'w') as f:
            f.write('#!/bin/bash\n')#module load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2=p
                f.write('$SCHRODINGER/run python {} {} {} {} {} {}\n'.format(script, l1, l2, st, glide_dir, out_dir))
            f.write('wait\n')
        os.system('sbatch --tasks=4 --cpus-per-task=1 --ntasks-per-socket=2 -p rondror -t 1:00:00 {}_in.sh'.format(i))

    os.chdir('../../../..')




