#!/bin/env run
import os
import click
import utils
import config

@click.group()
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/negative')
@click.option('--stats', default='/oak/stanford/groups/rondror/users/jpaggi/stats')
@click.pass_context
def main(ctx, data, stats):
	paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
			 'DATA': data, 'STATS': stats}
	paths.update(config.PATHS)
	paths = utils.resolve(paths)
	ctx.obj = paths

@main.command()
@click.argument('task')
@click.argument('stats_version', default='default')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def prepare(paths, task, stats_version, proteins):
	import dock.prepare_all
	params = config.STATS[stats_version]
	proteins = list(proteins)
	if not proteins:
		proteins = utils.get_proteins(paths, [])
	dock.prepare_all.main(params, paths, task, list(proteins))

# import ifp.fp
# import mcss.mcss
# import score.controller
# import score.scores
# import score.statistics
# elif mode == 'ifp':
#     ifp.fp.FP(args)
# elif mode == 'mcss':
#     mcss.mcss.main(args)
# elif mode == 'score_controller':
#     score.controller.main(args)
# elif mode == 'score':
#     score.scores.main(args)
# elif mode == 'statistics':
#     score.statistics.main(args)
# else:
#     print('Invalid arguments. Doing nothing.')


main()