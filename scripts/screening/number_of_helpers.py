from screening import *

threshes = np.linspace(0.2, 0.5, 13)

conditions = [(0,  np.linspace(1.00, 2.50, 7)),
			  (1,  np.linspace(1.00, 2.50, 7)),
			  (3,  np.linspace(0.75, 2.25, 7)),
			  (5,  np.linspace(0.75, 2.25, 7)),
			  (10, np.linspace(0.50, 2.0, 7)),
			  (15, np.linspace(0.50, 2.0, 7))]
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
perfs.to_csv('combind_search.csv')
