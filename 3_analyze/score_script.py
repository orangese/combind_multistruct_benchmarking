import os
import sys

#import numpy as np
#from multiprocessing import Pool

from containers import Dataset
from score_query import ScoreQuery

sys.path.append('../1_dock')
from parse_chembl import load_chembl

data = '/scratch/PI/rondror/jbelk'
chembl_dir = '{}/CHEMBL'.format(data)

dock_st = {'B1AR':'2VT4','AR_final':'3B67','TRPV1':'3J5Q'}
all_ligs = {p:[l.split('.')[0] for l in os.listdir('{}/{}/unique_ligands'.format(data, p)) 
               if l != st+'_lig.mae'] for p,st in dock_st.items()}

glide_dir = 'glide12'
cross = True
ifp_dir = 'ifp/ifp13'

w = {
    #1:2, # halogen bond
    2:1, # hydrogen bond
    3:1, # hydrogen bond
    4:1, # salt bridge
    6:0, # pipi
    #7:0, # picat
    #8:0, # picat
    10:0.01 # hydrophobic
}

#data_const = 0.5
num_poses = 25

prot = sys.argv[1]
query = sys.argv[2] # lig

#out_dir = 'scores3'
#all_lam = [i/10.0 for i in range(1,11)]
#all_w10 = [i/100.0 for i in range(6)]
#all_n = [i for i in range(21)]
#all_w56 = [(0,0),(0,1),(1,0)]

out_dir = 'scores5'
all_lam = [i/100.0 for i in range(5,16)]
all_w10 = [i/200.0 for i in range(4)]
all_n = [i for i in range(51)]
all_w56 = [(0,0),(0,1),(1,0)]

lam_i = int(sys.argv[3])
w10_i = int(sys.argv[4])
#n_i = int(sys.argv[5])
w56_i = int(sys.argv[5])

data_const = all_lam[lam_i]
w[10] = all_w10[w10_i]
#num_ligs = all_n[n]
w[5] = all_w56[w56_i][0]
w[6] = all_w56[w56_i][1]

d = Dataset([prot], data, dock_st, chembl_dir)
d.load_docking(glide_dir, ifp_dir, cross)
d.assign_weights(w)

lig_objs = d.all_proteins[prot].docking[(glide_dir, cross)].ligands
all_helper_ligs = [l for l in lig_objs if l[:6] == 'CHEMBL']
all_helper_ligs.sort(key = lambda x: lig_objs[x].ki)

description = 'best_affinity'
settings = 'prot={},dock_st={},glide_dir={},ifp_dir={},lambda={},num_poses={},w={}'.format(prot, dock_st[prot], glide_dir, ifp_dir, data_const, num_poses, w)

os.system('mkdir -p {}/{}/{}'.format(data, prot, out_dir))
os.system('mkdir -p {}/{}/{}/{}'.format(data, prot, out_dir, query))

def score_helper(num_ligs):
    f_name = '{}-{}-{}-{}.txt'.format(lam_i, w10_i, num_ligs, w56_i)
    f_path = '{}/{}/{}/{}/{}'.format(data, prot, out_dir, query, f_name)

    if os.path.exists(f_path): return#continue

    l_i = all_helper_ligs[:num_ligs]
    sq = ScoreQuery(query, l_i, lig_objs, data_const, num_poses, w)
    sq.score_all()

    out_ligs = [query] + sorted(l_i)

    with open(f_path, 'w') as f:
        f.write(description + '\n')
        f.write(settings + '\n')
        f.write(','.join(out_ligs) + '\n')
        for pi in range(sq.num_poses[query]):
            f.write(','.join([str(sq.pose_neighbors[pi][l].rank) for l in out_ligs]) + '\n')

#    return settings_tuple

#score_helper((query, lam_i, w10_i, w56_i))
#all_settings = [(query, lam_i, w10_i, w56_i, n) for n in all_n]

for n in all_n:
    score_helper(n)

#pool = Pool(int(os.environ.get("SLURM_NTASKS", 2)))
#print prot, all_ligs[prot]
#for q in pool.imap_unordered(score_helper, all_settings):
#    pass#print q, 'done'



