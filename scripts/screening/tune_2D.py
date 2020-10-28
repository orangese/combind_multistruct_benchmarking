from screening import *

threshes = [0.25, 0.275, 0.3]
weights_2D = np.linspace(10, 80, 15)
helpers = [0, 1, 3, 5, 10, 15]

perfs = []
for thresh in threshes:
	for helper in helpers:
		for weight in weights_2D:
		    print(thresh, helper, weight)
		    perf = plot_cat_all(thresh, thresh, weight_2D=weight, helpers=helper,
		                        helpers_sim=15, plot=False,
		                        scores=['COMBIND', 'GLIDE', 'COMBIND+2D', 'GLIDE+2D'])
		    perf['weight'] = weight
		    perf['helpers'] = helper
		    perf['thresh'] = thresh
		    perfs += [perf]

perfs = pd.concat(perfs).set_index(['thresh', 'helpers', 'weight', 'protein']).sort_index()
perfs.to_csv('tune_2D.csv')
