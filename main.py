#!/bin/sh
if "true" : '''\'
then
exec "$SCHRODINGER/run" "$0" "$@"
exit 127
fi
'''

import sys
# import dock.prepare_all
# import ifp.fp
# import mcss.mcss
import score.controller
import score.scores
import score.statistics
from settings import stats


mode = sys.argv[1]


if mode == 'prepare':
        dock.prepare_all.main(sys.argv[1:])

elif mode == 'ifp':
	ifp.fp.FP(sys.argv[1:])

elif mode == 'score_controller':
	score.controller.main(sys.argv[2:])

elif mode == 'mcss':
	mcss.mcss.main(sys.argv[1:])

elif mode == 'score':
	score.scores.main(sys.argv[1:])

elif mode == 'statistics':
	score.statistics.main(sys.argv[2:])

else:
	print('Invalid arguments. Doing nothing.')


