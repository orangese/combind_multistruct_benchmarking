import sys
import dock.prepare_all
import ifp.fp
import mcss.mcss
import score.controller
import score.scores
import score.statistics


if sys.argv[1] == 'prepare':
	dock.prepare_all.main(sys.argv[1:])
elif sys.argv[1] == 'ifp':
	ifp.fp.FP(sys.argv[1:])
elif sys.argv[1] == 'score_controller':
	score.controller.main(sys.argv[1:])
elif sys.argv[1] == 'mcss':
	mcss.mcss.main(sys.argv[1:])
elif sys.argv[1] == 'statistics':
	score.statistics.main(sys.argv[1:])
elif sys.argv[1] == 'score':
	score.scores.main(sys.argv[1:])

