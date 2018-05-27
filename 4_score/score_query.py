import os
import sys
import numpy

sys.path.append(os.getcwd())
from settings import *

code_path, q, s, p = sys.argv
code_path = '/'.join(code_path.split('/')[:-2])

sys.path.append(code_path+'/3_analyze')
sys.path.append(code_path+'/1_dock')

from containers import Dataset
from statistics import Statistics
from prob_opt import PredictStructs

fpath = '{}/{}-to-{}.sc'.format(os.getcwd(),q,s)

stats_data = Dataset(stats_prots, data_dir, glide_dir, ifp_dir, mcss_dir)
stats_data.load({p:prot.lm.get_docked(pdb_only=True) for p,prot in stats_data.proteins.items()})

stats = Statistics(stats_data, stats_prots, num_stats_ligs, num_stats_poses, features, smooth, normalize)

predict_data = Dataset([p], data_dir, glide_dir, ifp_dir, mcss_dir)
prot = predict_data.proteins[p]

chembl_ligs = prot.lm.get_similar(q, chembl_file, num=num_pred_chembl, mcss_sort=mcss_sort, struct=s)
predict_data.load({p:[q]+chembl_ligs},{p:[s]})

ps = PredictStructs(prot.docking[s], stats.evidence, features, num_poses, t)
   
best_cluster, all_scores, all_rmsds = ps.max_posterior([q]+chembl_ligs, restart=15, sampling=3)
result = ps.log_posterior(best_cluster)

us_top = best_cluster[q]#numpy.argmax(ps.score_query(q, best_cluster))

with open(fpath,'a') as f:
    f.write('{},{}\n'.format(q, us_top))
    for lig, pose in sorted(best_cluster.items()):
        if lig == q: continue
        f.write('{},{}\n'.format(lig, pose))

    f.write('max_score,{}\n'.format(result))


