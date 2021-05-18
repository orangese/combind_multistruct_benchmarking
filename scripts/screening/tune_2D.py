"""
sbatch -p rondror --wrap="python scripts/screening/tune_2D.py ef1"  -t 6:00:00
sbatch -p rondror --wrap="python scripts/screening/tune_2D.py logroc"  -t 6:00:00
"""

from screening import *
import sys

metric = sys.argv[1]

threshes = [('0-2', 0, 0.2), ('2-3', 0.2, 0.3)]
helpers = [0, 1, 5, 10]

scores =  ['COMBIND', 'GLIDE',
		   'COMBIND_GLIDE', 'GLIDE_COMBIND',
	       '2D_mean', 'SHAPE_mean',
	       'Z:COMBIND+2D+SHAPE', 'Z:GLIDE+2D+SHAPE',
	       'N:COMBIND+2D+SHAPE', 'N:GLIDE+2D+SHAPE',
	       'R:COMBIND+2D+SHAPE', 'R:GLIDE+2D+SHAPE',
	       'N:2D+SHAPE'
	       ]

perfs = []
for (thresh, lower, upper) in threshes:
	for helper in helpers:
	    print(thresh, helper)
	    perf = plot_cat_all(upper, upper, active_min_cut=lower, helpers=helper,
	                        helpers_sim=10, scores=scores, metric=metric,
	                        minimum_active=10)
	    perf['helpers'] = helper
	    perf['thresh'] = thresh
	    perfs += [perf]

perfs = pd.concat(perfs).set_index(['thresh', 'helpers', 'protein']).sort_index()
perfs.to_csv('{}.csv'.format(metric))
