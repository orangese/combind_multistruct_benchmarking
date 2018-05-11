import os
import sys
import numpy

sys.path.append(os.getcwd())
from settings import *

code_path, q, p = sys.argv
code_path = '/'.join(code_path.split('/')[:-2])

sys.path.append(code_path+'/3_analyze')
sys.path.append(code_path+'/1_dock')

from containers import Dataset
from statistics import Statistics
from prob_opt import PredictStructs

stats_data = Dataset(stats_prots, data_dir, glide_dir, ifp_dir, mcss_dir)
stats_data.load({p:prot.lm.pdb for p,prot in stats_data.proteins.items()})

stats = Statistics(stats_data, stats_prots, num_stats_ligs, num_stats_poses, features, smooth)

predict_data = Dataset([p], data_dir, glide_dir, ifp_dir, mcss_dir)
predict_data.load({p:prot.lm.pdb for p, prot in predict_data.proteins.items()})

prot = predict_data.proteins[p]
ps = PredictStructs(prot.docking[prot.lm.default_st], stats.evidence, features, num_poses, t)
    
chembl_ligs = prot.lm.get_similar(q, chembl_file, num=num_pred_chembl, mcss_sort=mcss_sort)
predict_data.load({p:[q]+chembl_ligs})
        
best_cluster, en_landscape = ps.max_posterior(chembl_ligs)
result = ps.joint_posterior(best_cluster)

us_top = numpy.argmax(ps.score_query(q, best_cluster))

print '{},{}'.format(q, us_top)
for lig, pose in sorted(best_cluster.items()):
    print '{},{}'.format(lig, pose)

print 'max_score,{}'.format(result[0])


