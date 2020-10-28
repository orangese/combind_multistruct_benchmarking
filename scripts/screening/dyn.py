from screening import *

threshes = np.linspace(0.2, 0.5, 13)

conditions = [(0,  [1.0]),
			  (1,  [1.0]),
			  (3,  [1.0]),
			  (5,  [1.0]),
			  (10, [1.0]),
			  (15, [1.0])]
perfs = []
for thresh in threshes:
	for helpers, weights in conditions:
		for weight in weights:
		    print(thresh, helpers, weight)
		    perf = plot_cat_all(thresh, thresh, weight_combind=weight, helpers=helpers,
		                        helpers_sim=15, plot=False,
		                        scores=['COMBIND', 'GLIDE', 'COMBIND_GLIDE',
		                                'GLIDE_COMBIND', '2D_mean', 'SHAPE_mean'])
		    perf['weight'] = weight
		    perf['helpers'] = helpers
		    perf['thresh'] = thresh
		    perfs += [perf]

perfs = pd.concat(perfs).set_index(['thresh', 'helpers', 'weight', 'protein']).sort_index()
perfs.to_csv('combind_dyn.csv')
