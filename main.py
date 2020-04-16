#!/bin/env run
import os
import click
import utils
import config

@click.group()
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/negative')
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

@main.command()
@click.option('--plot', is_flag=True)
@click.argument('stats_root')
@click.argument('struct')
@click.argument('protein')
@click.argument('queries', nargs=-1)
@click.pass_obj
def score(paths, stats_root, struct, protein, queries, plot):
	import score.scores
	queries = list(queries)
	score.scores.main(paths, config.FEATURE_DEFS,
	                  stats_root, struct, protein, queries, plot)

@main.command()
@click.option('--merged_root')
@click.option('--plot')
@click.argument('stats_version')
@click.argument('stats_root')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def statistics(paths, stats_version, stats_root, proteins, merged_root, plot):
	import score.statistics
	params = config.STATS[stats_version]
	proteins = list(proteins)
	score.statistics.main(params, paths, config.FEATURE_DEFS, stats_root,
	                      proteins, merged_root, plot)

@main.command()
@click.argument('ifp_version')
@click.argument('input_file')
@click.argument('output_file')
@click.argument('poses', type=int)
def ifp(ifp_version, input_file, output_file, poses):
	import ifp.ifp
	ifp.ifp.IFP(config.IFP[ifp_version], input_file, output_file, poses)

@main.command()
@click.argument('mode')
@click.argument('ligand1')
@click.argument('ligand2')
@click.argument('ligand1_path')
@click.argument('ligand2_path')
@click.argument('init_file')
@click.argument('mcss_types_file')
@click.argument('rmsd_file', default='')
@click.argument('max_poses', type=int, default=0)
@click.argument('mcss_string_rep', default='')
def mcss(mode, ligand1, ligand2, ligand1_path, ligand2_path,
         init_file, mcss_types_file, rmsd_file, max_poses,
         mcss_string_rep):
	import mcss.mcss
	mcss.mcss.main(mode, ligand1, ligand2, ligand1_path, ligand2_path,
                   init_file, mcss_types_file, rmsd_file, max_poses,
                   mcss_string_rep)

main()
