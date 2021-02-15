from screening import *

threshes = np.linspace(.20, 0.30, 3)
helpers = [0, 1, 3, 5, 10]

perfs = []
for thresh in threshes:
	thresh = round(thresh, 2)
	for helper in helpers:
	    print(thresh, helper)
	    perf = plot_cat_all(thresh, thresh, helpers=helper,
	                        helpers_sim=10, plot=False,
	                        scores=['COMBIND', 'GLIDE', '2D_mean', 'SHAPE_mean',
	                        		'Z:COMBIND+2D', 'Z:GLIDE+2D',
	                        		'Z:COMBIND+SHAPE', 'Z:GLIDE+SHAPE',
	                        		'Z:COMBIND+2D+SHAPE', 'Z:GLIDE+2D+SHAPE',
	                        		'N:COMBIND+2D', 'N:GLIDE+2D',
	                        		'N:COMBIND+SHAPE', 'N:GLIDE+SHAPE',
	                        		'N:COMBIND+2D+SHAPE', 'N:GLIDE+2D+SHAPE',
	                        		'N:2D+SHAPE', 'Z:2D+SHAPE'
	                        		])
	    perf['helpers'] = helper
	    perf['thresh'] = thresh
	    perfs += [perf]

perfs = pd.concat(perfs).set_index(['thresh', 'helpers', 'protein']).sort_index()
perfs.to_csv('tune_2D_EF.csv')
