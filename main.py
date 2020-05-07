#!/bin/env run
import os
import click
import utils
import config

@click.group()
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/negative')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
@click.pass_context
def main(ctx, data, ligands):
    paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
             'DATA': data,
             'PDB': ligands}
    paths.update(config.PATHS)
    paths = utils.resolve(paths)
    ctx.obj = paths

@main.command()
@click.option('--stats_version', default='default')
@click.option('--struct', default=None)
@click.argument('task')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def prepare(paths, task, stats_version, proteins, struct):
    """
    Prep structures and ligands; docking; and featurization.
    """
    import dock.prepare_all
    params = config.STATS[stats_version]
    proteins = list(proteins)
    if not proteins:
        proteins = utils.get_proteins(paths, [])
    dock.prepare_all.main(params, paths, task, list(proteins), struct)

@main.command()
@click.option('--plot', is_flag=True)
@click.option('--stats_root', default='{}/stats_data/rd1'.format(os.environ['COMBINDHOME']))
@click.option('--struct', default=None)
@click.option('--fname', default='pdb.sc')
@click.option('--xtal', multiple=True)
@click.argument('protein')
@click.argument('queries', nargs=-1)
@click.pass_obj
def score(paths, stats_root, struct, protein, queries, plot, fname, xtal):
    """
    Run ComBind!
    """
    import score.scores
    queries = list(queries)
    xtal = list(xtal)
    score.scores.main(paths, config.FEATURE_DEFS, stats_root, protein, queries,
                      xtal=xtal, fname=fname, struct=struct, plot=plot)

@main.command()
@click.option('--stats_root', default='{}/stats_data/rd1'.format(os.environ['COMBINDHOME']))
@click.option('--struct', default=None)
@click.argument('protein')
@click.argument('cluster')
@click.argument('queries', nargs=-1)
@click.pass_obj
def screen(paths, stats_root, struct, protein, queries, cluster):
    """
    Run ComBind virtual screening!
    """
    import score.scores
    queries = list(queries)
    score.scores.screen(paths, config.FEATURE_DEFS,
                        stats_root, struct, protein, queries, cluster)

@main.command()
@click.option('--merged_root')
@click.option('--stats_version', default='default')
@click.argument('stats_root')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def statistics(paths, stats_version, stats_root, proteins, merged_root):
    """
    Fit statistics to a set of protein–ligand complexes.
    """
    import score.statistics
    params = config.STATS[stats_version]
    proteins = list(proteins)
    if not proteins:
        proteins = utils.get_proteins(paths, [])
    score.statistics.compute(params, paths, config.FEATURE_DEFS, stats_root,
                             proteins, merged_root)

@main.command()
@click.argument('merged')
def plot_statistics(merged):
    import score.statistics
    score.statistics.plot(merged)

@main.command()
@click.argument('ifp_version')
@click.argument('input_file')
@click.argument('output_file')
@click.argument('poses', type=int)
def ifp(ifp_version, input_file, output_file, poses):
    """
    Compute interaction fingerprints. (For internal use.)
    """
    import ifp.ifp
    ifp.ifp.ifp(config.IFP[ifp_version], input_file, output_file, poses)

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
    """
    Compute substructure similiarity. (For internal use.)
    """
    import mcss.mcss
    mcss.mcss.main(mode, ligand1, ligand2, ligand1_path, ligand2_path,
                   init_file, mcss_types_file, rmsd_file, max_poses,
                   mcss_string_rep)

main()
