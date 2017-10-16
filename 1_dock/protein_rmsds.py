import os
from multiprocessing import Pool
script = '/share/PI/rondror/jbelk/combind/1_dock/rmsd_by_residue.py'

schro = os.environ.get("SCHRODINGER", None)

def calc_rmsds():
    os.system('rm -rf rmsds')
    os.system('mkdir -p rmsds')

    all_prot = sorted([i.split('_')[0] for i in os.listdir('final_ligands')])
    all_pairs = []
    
    empty = 0
    for f in os.listdir('rmsds'):
        if f.split('.')[-1] == 'csv' and os.path.getsize('rmsds/{}'.format(f)) == 0:
            empty += 1
            os.system('rm -f rmsds/{}'.format(f))
    print 'deleted {} empty rmsd files'.format(empty)

    ref = all_prot[0]
    for i in all_prot[1:]:
        #for j in range(i+1, len(all_prot)):
        if '{}-{}.csv'.format(ref, i) not in os.listdir('rmsds'):
            all_pairs.append((ref, i))
    print len(all_pairs)
    if len(all_pairs) == 0: return
    
    pool = Pool(int(os.environ.get("SLURM_NTASKS", 5)))
    
    for s in pool.imap_unordered(helper, all_pairs):
        print s

def helper(s):
    with open('rmsds/{}-{}.csv'.format(s[0],s[1]), 'w+') as f:
        pass
    os.system('{}/run {} -o rmsds/{}-{}.csv final_proteins/{}.mae final_proteins/{}.mae -d 5'.format(schro, script, s[0], s[1], s[0], s[1]))
