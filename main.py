#!/bin/env run
import sys
import dock.prepare_all
import ifp.fp
import mcss.mcss
import score.controller
import score.scores
import score.statistics
from settings import stats

mode, args = sys.argv[1], sys.argv[2:]
if mode == 'prepare':
    dock.prepare_all.main(args)
elif mode == 'ifp':
	ifp.fp.FP(args)
elif mode == 'mcss':
	mcss.mcss.main(args)
elif mode == 'score_controller':
	score.controller.main(args)
elif mode == 'score':
	score.scores.main(args)
elif mode == 'statistics':
	score.statistics.main(args)
else:
	print('Invalid arguments. Doing nothing.')
